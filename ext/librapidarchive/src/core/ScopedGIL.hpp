#pragma once

#include <exception>
#include <optional>
#include <stdexcept>
#include <vector>

#include <Python.h>


namespace rapidgzip
{
[[nodiscard]] bool
pythonIsFinalizing()
{
    #if ( PY_MAJOR_VERSION != 3 ) || ( PY_MINOR_VERSION < 8 )
        return false;
    #elif PY_MINOR_VERSION < 13
        return _Py_IsFinalizing();
    #else
        return Py_IsFinalizing();
    #endif
}


class PythonExceptionThrownBySignal :
    public std::runtime_error
{
public:
    PythonExceptionThrownBySignal() :
        std::runtime_error( "An exception has been thrown while checking the Python signal handler." )
    {}
};


class ScopedGIL
{
public:
    struct GILState
    {
        bool locked;
        bool exists;
    };

public:
    explicit
    ScopedGIL( bool doLock )
    {
        m_referenceCounters.emplace_back( apply( { /* .locked = */ doLock, /* .exists = */ true } ) );
    }

    ~ScopedGIL() noexcept
    {
        if ( m_referenceCounters.empty() ) {
            std::cerr << "Logic error: It seems there were more unlocks than locks!\n";
            std::terminate();
        }

        apply( m_referenceCounters.back() );
        m_referenceCounters.pop_back();
    }

    ScopedGIL( const ScopedGIL& ) = delete;
    ScopedGIL( ScopedGIL&& ) = delete;
    ScopedGIL& operator=( const ScopedGIL& ) = delete;
    ScopedGIL& operator=( ScopedGIL&& ) = delete;

private:
    /**
     * @return the old GIL state.
     */
    GILState
    apply( const GILState targetState ) noexcept
    {
        const auto doLock = targetState.locked;
        if ( !doLock && pythonIsFinalizing() ) {
            /* No need to unlock the GIL if it doesn't exist anymore. Should not matter what we return here. */
            return { /* .locked = */ false, /* .exists = */ false };
        }

        if ( targetState.locked && !targetState.exists ) {
            std::cerr << "Invalid GIL target state, which should be locked but not exist at the same time!\n";
            std::terminate();
        }

        /**
         * I would have liked a GILMutex class that can be declared as a static thread_local member but
         * on Windows, these members are initialized too soon, i.e., at static initialization time instead of
         * on first usage, which leads to bugs because PyGILState_Check will return 0 at this point.
         * Therefore, use block-scoped thread_local variables, which are initialized on first pass as per the standard.
         * @see https://stackoverflow.com/a/49821006/2191065
         */
        static thread_local bool isLocked{ PyGILState_Check() == 1 };

        /** Used for locking non-Python threads. */
        static thread_local std::optional<PyGILState_STATE> lockState{};
        /** Used for unlocking and relocking the Python main thread. */
        static thread_local PyThreadState* unlockState{ nullptr };

        /* When Python is finalizing, we might get our acquired GIL rugpulled from us, meaning isLocked=true but
         * PyGILState_Check() is 0 / unlocked!
         * Python 3.10 has _Py_IsFinalizing, 3.13 has Py_IsFinalizing, however on Python 3.6, these are missing. */
        if ( pythonIsFinalizing() || ( isLocked && ( PyGILState_Check() == 0 ) ) ) {
            if ( ( PyGILState_Check() == 1 ) && lockState.has_value() ) {
                PyGILState_Release( *lockState );
                lockState.reset();
            }
            std::cerr << "Detected Python finalization from running rapidgzip thread.\n"
                         "To avoid this exception you should close all RapidgzipFile objects correctly,\n"
                         "or better, use the with-statement if possible to automatically close it.\n";
            std::terminate();
        }

        const auto wasLocked = isLocked;
        if ( isLocked == doLock ) {
            return { /* .locked = */ wasLocked, /* .exists = */ true };
        }

        /* PyPy did not have PyGILState_GetThisThreadState, so we have to use a workaround.
         * https://github.com/pypy/pypy/issues/5302 */
        PyThreadState* threadState{ nullptr };
        #ifdef PYPY_VERSION_NUM
            threadState = _PyThreadState_UncheckedGet();
        #else
            threadState = PyGILState_GetThisThreadState();
        #endif

        const auto gilExists = threadState != nullptr;
        if ( doLock ) {
            if ( gilExists ) {
                /* Calling PyEval_RestoreThread with PyGILState_GetThisThreadState is what PyGILState_Ensure does
                 * internally, so it should be fine to have that as fallback. */
                PyEval_RestoreThread( unlockState == nullptr ? threadState : unlockState );
                unlockState = nullptr;
            } else {
                /* This should only happen on spawned C-threads when they try to call into Python the very first
                 * time. All recursive calls into Python should avoid PyGILState_Ensure/_Release because destroying
                 * and recreating the thread state can lead to use-after-free bugs!
                 * Therefore, this should only happen when m_referenceCounters is empty.
                 * The value in lockState should always be PyGILState_UNLOCKED. */
                lockState.emplace( PyGILState_Ensure() );
            }
        } else {
            if ( !targetState.exists && lockState.has_value() ) {
                /* https://github.com/python/cpython/blob/5334732f9c8a44722e4b339f4bb837b5b0226991/Python/pystate.c
                 *   #L2833
                 * PyGILState_Release basically does this: if (oldstate == PyGILState_UNLOCKED) PyEval_SaveThread();
                 * Furthermore, it DESTRUCTS the thread state when the counter reaches 0! This can be bad because
                 * it invalidates all thread state pointers, which may lead to bugs if some calling function still
                 * has such a pointer and tries to use it after we return! */
                PyGILState_Release( *lockState );
                lockState.reset();
            } else {
                unlockState = PyEval_SaveThread();
            }
        }

        isLocked = doLock;
        return { /* .locked = */ wasLocked, /* .exists = */ gilExists };
    }

private:
    inline static thread_local std::vector<GILState> m_referenceCounters;
};


class ScopedGILLock :
    public ScopedGIL
{
public:
    ScopedGILLock() :
        ScopedGIL( true )
    {}
};


class ScopedGILUnlock :
    public ScopedGIL
{
public:
    ScopedGILUnlock() :
        ScopedGIL( false )
    {}
};


void
checkPythonSignalHandlers()
{
    const ScopedGILLock gilLock;

    /**
     * @see https://docs.python.org/3/c-api/exceptions.html#signal-handling
     * > The function attempts to handle all pending signals, and then returns 0.
     * > However, if a Python signal handler raises an exception, the error indicator is set and the function
     * > returns -1 immediately (such that other pending signals may not have been handled yet:
     * > they will be on the next PyErr_CheckSignals() invocation).
     */
    for ( auto result = PyErr_CheckSignals(); result != 0; result = PyErr_CheckSignals() ) {
        if ( PyErr_Occurred() != nullptr ) {
            throw PythonExceptionThrownBySignal();
        }
    }
}
}  // namespace rapidgzip
