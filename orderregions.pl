#!/usr/bin/perl -w
#
# orderregions.pl
# This script takes the collated csv output from the fragcollater and the raw
# blast output and orders the regions as they would appear in the protein. It
# also provides the position in the target protein.
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

# list of organisms of interest
my @organism = qw( micro volvox auxeno cocco ostreo_t ostreo_l chlamy chlorella selag physco );


foreach my $organism (@organism) {
    
    print "[$organism] Processing ";
    # get list of csv files
    my @fragblasts = <$organism/fragblast/*.csv>;
    my @blasts = qw( BLASTP DELTA PSIBLAST PSSM );
    
    foreach my $fragblast (@fragblasts) {
        
        next if $fragblast =~ /FULL/;
        
        my ($file) = $fragblast =~ /$organism\/fragblast\/(.+)\.csv/;
        my @lines = read_file($fragblast);
        #retreive header from csv
        my $header = shift @lines;
        chomp $header;
        # order header regions
        my @regions = split(/,/,$header);
        my @sorted = orderRegions(@regions);
        my %globalhits;
        
        foreach my $region (@sorted) {
            
            #skip if the extra field
            next if $region eq '';
            # get associated blast outputs for the protein region
            my $rootfilename = $file."___".$region."____";
            my @blastfiles = <$organism/blastout/$rootfilename*.txt>;
            
            # grabs hit info for each hit in a region
            foreach my $blastout (@blastfiles) {
                
                my ($basefile,$blasttype) = $blastout =~ /$organism\/blastout\/(.+)____(.+).txt/;
                my $filecontent = read_file($blastout);
                #gets the position of the hits in the target organism
                %{$globalhits{$basefile}{$blasttype}} = processFile($filecontent);
            }
            
        }
        
        my @output;
        open ORDERED, ">", "$organism/ordered/".$file.".ordered";
        print ORDERED ",",join(",",@sorted),"\n";
        
        # builds output file rows
        foreach my $row (@lines) {
            
            # grabs protein version id from csv file row puts it on the front of the output row.
            chomp $row;
            my @row = split(/,/,$row);
            my @outrow;
            push(@outrow,$row[0]);
            
            my $length;
            my $hitsomething = 0;
            
            # cycles through query protein regions and blast types to collate the hits for each region.
            for (my $i = 0; $i < scalar @sorted; $i++) {
                
                my $key = $file."___".$sorted[$i];
                my @regionhits;
                # collates blast hits for each reagion from the globalhits table.
                foreach my $blasttype (@blasts) {
                    
                    if (exists $globalhits{$key}{$blasttype}{$row[0]} ){
                        $hitsomething = 1;
                        $length = $globalhits{$key}{$blasttype}{$row[0]}{'GENELENGTH'}; 
                        #my $hit = $globalhits{$key}{$blasttype}{$row[0]};
                        push(@regionhits,$blasttype.':'.$globalhits{$key}{$blasttype}{$row[0]}{'EXPECT'}.':'.$globalhits{$key}{$blasttype}{$row[0]}{'START'}.'..'.$globalhits{$key}{$blasttype}{$row[0]}{'END'});
                    }
                    
                }
                
                push(@outrow,join("|",@regionhits));
                
            }
            
            # add protein length after the version id.
            splice @outrow, 1, 0, $length;
            #print the row only if there is a row so perl doesn't complain
            print ORDERED join(",",@outrow),"\n" if $hitsomething;
            push (@output,\@outrow);
            
        }
        
        
        print ".";
        
    }
    print "\n";
}
#sorts regions by residue positon
sub orderRegions {
    
    my %regions;
    shift @_;
    foreach (@_) {
        
        unless (/FULL/) {
            
            chomp;
            my ($start,$finish,$name) = /(\d+)\.\.(\d+)_(.+)/;
            die "$start already exists: $name" if exists $regions{$start};
            $regions{$start} = $_;
            
        }
        
    }
    
    return map {$regions{$_}} sort{ $a <=> $b } keys %regions;
    
}
# grabs hit data from raw blast output.
sub processFile {
    
    my $filecontent = shift;
    my %hits;
    foreach ($filecontent =~ m/(>[^>]+)/mg) {
        # match data points of interest from the file.
        my ($version, $expect) = $_ =~ /^>\s([^\s]+).+Expect\s=\s([^\,]+)/sm;
        my ($length) = $_ =~ m/^Length=(\d+)/ms; 
        my ($first) = split(/\n\n\n/);
        my ($start,$end) = $first =~ m/Sbjct\s\s(\d+)[^>]+\d*\s+[^\s]+\s+(\d+)$/mg;
        
        # stuff the data into a table
        $hits{$version}{'START'} = $start;
        $hits{$version}{'END'} = $end;
        $hits{$version}{'EXPECT'} = $expect;
        $hits{$version}{'GENELENGTH'} = $length;
        
    }
    return %hits;
}
#
sub psiblastSegmenter {
    my $file = shift;
    my @fragments = split(/Results from round \d\d?/,$file);
    return pop @fragments;
}

