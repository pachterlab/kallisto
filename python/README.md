The following are some notes that Nick Bray provided in an email describing the
code

----

- The way make\_equiv\_classes works is by iteratively building the equivalence
  classes as you add transcripts. For a given k-mer in your new transcript, if
  you've never seen that k-mer before, it belongs to a new equivalence class of
  k-mers that map just to that transcript. If you have seen it before, then
  whatever equivalence class it was in before is split into those k-mers that
  map to this new transcript and those that don't. Probably clumsy but...it got
  the job done. Of course the details of how this is done is probably
  completely uninteresting since I'm sure it's getting reimplemented anyway but
  one property of this is used later on: the equivalence classes are given
  consecutive id numbers as they're created and the order in which they're
  created here ensures that the ordering of those id numbers agrees with the
  lexicographic ordering of the count vectors. That means that if you take the
  minimum of two count vectors, it will be <= the equivalence class with the
  smaller id number. Now that I'm writing this, I'm asking myself "Does this
  actually help the program in any way?" and I'm not 100% sure. I thought at
  the time that it would make things at least a touch more efficient though.

- You'll notice that currently the equivalence classes and such are being
  recomputed every time you run a dataset rather than once per transcriptome.
  Pickling the results seemed to be unbelievably slow so I just gave up on that
  (especially since right now it's a pretty small fraction of the time taken).

- The equivalence classes give the counts of the number of times a k-mer occurs
  in each transcript but you also need to track the strand (well, maybe) so the
  way this works is equiv\_classes[i][2*t + s] is the number of times a k-mer in
  equivalence class i occurs in transcript t on strand s (1 if reverse,
  0 otherwise).

I say "well, maybe" because the fact that we're taking the min here means that
we could probably just ignore strand and everything will work out fine. It
would probably speed things up a bit too, so that might be something to try at
some point.

- Oh, so one thing is that right now this is written only for paired end reads.
  It shouldn't be a big deal to change that but it does seem that as far as
  getting TPM right, having accurate effective lengths is pretty important and
  I don't have any good idea for how you'd do that with single end reads.
  (Currently there's a script for generating a fragment length distribution
  based on a smallish number of read alignments. I'll be sending that along
  later.)

- calc_trans_ec_weights basically gives you something like P(read|transcript)
  for a given transcript and a read in a given equivalence class. Right now
      this is just being done by (count in the equivalence class
      vector)/(effective length of transcript). There's a fancier way of doing
      this (accessed by -st and currently *very very* slow) which I was really
      anticipating would be better but that doesn't seem to be the case at the
      moment, unfortunately.

- process_reads actually takes the reads and assigns them to equivalence
  classes. Its operation is actually not that complicated once you hide some
  ugly stuff inside the memoized_min function.

- calc_alpha actually does the EM algorithm. It works by splitting each
  equivalence class between the transcripts its reads could have come from (we
  can definitely ignore strand at this point) and I think the only particular
  spin on the standard RNA-seq EM is that it calculates the denominator for
  each splitting first. That's probably...completely opaque. But maybe we
  should just talk in person tomorrow anyway.
