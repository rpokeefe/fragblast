#!/usr/bin/perl -w
#
# fragcollater.pl
# This script parses the raw blast output files and organizes the output by
# gene hit version number into csv files. Columns are the query regions, rows
# are the protein hit version number.
#
# This script and the accompanying scripts were written for one time use
# to investigate a specific question, but if you want to branch and develop
# further GPL license info follows.
#                             ****
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#############################################################################

use strict;
use File::Slurp;

#list of investigated organisms
my @organism = qw( micro volvox auxeno cocco ostreo_t ostreo_l chlamy chlorella selag physco );

foreach my $organism (@organism) {
    
    # get list of blast output files
    my @files = <$organism/blastout/*.txt>;
    my %file;
    
    # Loop parses and extracts hits and evalues from raw blast output.
    foreach  (@files) {
        
        print "[$organism] Processing: $_\n";
        my $blastout = read_file($_);
        
        #parse filename for context
        my ($fullprotein, $region, $action) = $_ =~ m/$organism\/blastout\/(.+___.+)___(.+)____(.+)\.txt/;
        if ($action =~ /PSIBLAST/) {
            
            #psiblast have a segment per iteration and we only want the last iteration.
            $blastout = psiblastSegmenter($blastout);
        }
        #grab lines that match result summary lines.
        my @matches = $blastout =~ m/^\s\s([^\s]{10}).+\s{2,}\d+\.?\d+\s{2,}([\d\.e\-]+)\s*$/gm;
        
        #grab each id and evalue from the matches and dumpt them into our datastructure.
        while (my ($id, $expected) = splice(@matches,0,2)) {
            
            #this creates a counter for the hit protein for sorting later.
            unless (defined $file{$fullprotein}{'REGION'}{$region}{$id}{'EXPECT'}) {
                $file{$fullprotein}{'COUNT'}{$id}++; 
            }
            # initialize the region with an evalue if it doesn't exist or update
            # the evalue if the new value has a lower value.
            unless (defined $file{$fullprotein}{'REGION'}{$region}{$id}{'EXPECT'}) {
                $file{$fullprotein}{'REGION'}{$region}{$id}{'EXPECT'} = $expected;
            } elsif ($file{$fullprotein}{'REGION'}{$region}{$id}{'EXPECT'} > $expected) {
                $file{$fullprotein}{'REGION'}{$region}{$id}{'EXPECT'} = $expected;
            }
            
        }
        
    }
    
    # loop generates an output file for each query protein's blast results for a specific organism.
    foreach my $fullprotein (keys %file) {
        
        # get the regions of the query protein
        my @regions = keys %{$file{$fullprotein}{'REGION'}};
        # add a space to the front of the list to shift the columns in the file header.
        push(@{$file{$fullprotein}{'FILE'}[0]},'',@regions);
        
        # sort the hit proteins by the number of regions where hits occured. Loop constructs the data grid that will be dumped into the file.
        foreach my $id (sort {$file{$fullprotein}{'COUNT'}{$b} <=> $file{$fullprotein}{'COUNT'}{$a}} keys %{$file{$fullprotein}{'COUNT'}}) {
            
            #get the current row of the file matrix
            my $rowindex = scalar @{$file{$fullprotein}{'FILE'}};
            #set protein version id as first element of row
            $file{$fullprotein}{'FILE'}[$rowindex][0] = $id;
            
            
            for (my $i = 0; $i < scalar @regions; $i++) {
                
                # populate my file grid with available evalues or add a blank field.
                if (defined $file{$fullprotein}{'REGION'}{$regions[$i]}{$id}{'EXPECT'}) {
                    $file{$fullprotein}{'FILE'}[$rowindex][$i + 1] = $file{$fullprotein}{'REGION'}{$regions[$i]}{$id}{'EXPECT'};
                } else {
                    $file{$fullprotein}{'FILE'}[$rowindex][$i + 1] = '';
                }
                
            }
            
        }
        # create and open output file
        my $outfile = "$organism/fragblast/".$fullprotein.'.csv';
        print "[$organism] Writing $outfile \n";
        open OUTFILE, ">", $outfile or die $!;
        # dump matrix to file
        print OUTFILE join("\n",map { join(',',@{$_})} @{$file{$fullprotein}{'FILE'}});
        close OUTFILE;
        
    }
}

#grabs the summary from the last iteration of a psiblast.
sub psiblastSegmenter {
    my $file = shift;
    my @fragments = split(/Results from round \d\d?/,$file);
    return pop @fragments;
}




