#!/usr/bin/env perl

# Copyright (c) 2025 GilbertHan (gilberthan1011@gmail.com)
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.
#
# ALSO, IT WOULD BE NICE IF YOU LET ME KNOW YOU USED IT.

# This script is adapted to output TSS (Transcription Start Site)
# in BED6 format: chrom, start, end, transcript_id, 0, strand

use Getopt::Long;

$in = shift @ARGV;

my $in_cmd =($in =~ /\.gz$/ ? "gunzip -c $in|" : $in =~ /\.zip$/ ? "unzip -p $in|" : "$in") || die "Can't open $in: $!\n";
open IN, $in_cmd;

my %tss_seen;
while (<IN>) {
    $gff = 2 if /^##gff-version 2/;
    $gff = 3 if /^##gff-version 3/;
    next if /^#/ && $gff;

    s/\s+$//;
    # 0-chr 1-src 2-feat 3-beg 4-end 5-scor 6-dir 7-fram 8-attr
    my @f = split /\t/;
    next unless @f >= 9;

    my $id;
    if ($gff) {
        ($id) = $f[8] =~ /\bID="([^"]+)"/;
        ($id) = $f[8] =~ /\bName=([^";]+)/ if !$id && $gff == 3;
    } else {
        ($id) = $f[8] =~ /transcript_id "([^"]+)"/;
    }
    next unless $id && $f[0];

    # only process feature types that define a transcript start
    # Use 'transcript' or 'mRNA' or 'miRNA'
    if ($f[2] eq 'transcript' || $f[2] eq 'mRNA' || $f[2] eq 'miRNA') {
        my ($chr, $beg, $end, $strand) = ($f[0], $f[3], $f[4], $f[6]);
        # BED is 0-based, half-open.
        my ($tss_start, $tss_end);

        if ($strand eq '+') {
            $tss_start = $beg - 1;
            $tss_end = $tss_start + 1;
        } else {
            $tss_end = $end;     # end is 1-based inclusive
            $tss_start = $end - 1;
        }

        # Remove duplicate TSS for same transcript if present
        next if $tss_seen{"$chr:$tss_start:$strand:$id"}++;
        print "$chr\t$tss_start\t$tss_end\t$id\t0\t$strand\n";
    }
}

close IN;
