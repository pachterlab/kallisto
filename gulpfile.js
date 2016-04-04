var gulp = require('gulp');
var exec = require('child_process').exec;

var cMakeCommand = 'cd build; cmake ..;';
var buildCommand = 'cd build; make;';

var indexCommand = 'build/src/kallisto' +
  ' index -i test/transcripts.kidx' +
  ' test/transcripts.fasta.gz';

var pairedEndCommand = 'build/src/kallisto' +
  ' quant -i test/transcripts.kidx' +
  ' -b 10' +
  ' -t 2' +
  ' -o test/paired_end' +
  ' test/reads_1.fastq.gz test/reads_2.fastq.gz';

var singleEndCommand = 'build/src/kallisto' +
  ' quant -i test/transcripts.kidx' +
  ' -b 10' +
  ' -t 2' +
  ' -l 200 -s 3' +
  ' -o test/single_end' +
  ' --single' +
  ' test/reads_1.fastq.gz';

console.log('build command: ' + buildCommand);

gulp.task('watch', function() {
  gulp.watch('src/*.cpp', ['build']);
  gulp.watch('src/*.h', ['build']);
  gulp.watch('src/*.hpp', ['build']);
});

gulp.task('build', ['watch'], function() {
  exec(buildCommand, function(error, standardOutput, standardError) {
    if (error) {
      console.error('There was an error: ' + error);
    }
    console.log(standardOutput);
    console.log(standardError);
  });
});

gulp.task('pairedEnd', ['build'], function() {
  exec(pairedEndCommand, function(error, standardOut, standardError) {
    if (error) {
      console.error('There was a pairedEnd error');
    }
    console.log(standardOut);
    console.log(standardError);
  });
});

gulp.task('singleEnd', ['build'], function() {
  exec(singleEndCommand, function(error, standardOut, standardError) {
    if (error) {
      console.error('There was a singleEnd error');
    }
    console.log(standardOut);
    console.log(standardError);
  });
});

// gulp.task('default', ['install', 'watch'], function() {});
// gulp.task('compile', ['build', 'watch'], function() {});
// gulp.task('pairedEnd', ['compile'], function() {});
// gulp.task('singleEnd', ['compile'], function() {});
gulp.task('default', ['pairedEnd', 'singleEnd'], function() {});
