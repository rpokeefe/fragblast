package GenBankEasy;
#############################################################################
# GenBankEasy.pm
# This package provides a simple parser for Genbank files to accommodate
# the fragmented blast methods used in this process.
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
use warnings;

sub new {
    
    my $class = shift;
    my $self = bless {}, $class;
    return $self;
}


sub loadRec {
    
    my $self = shift;
    my $rec = shift;
    
    my @rec = split(/\n/, $rec);
    my (%GB,$tagname);
    
    foreach my $line (@rec) {
        
        if ($line =~ /^([A-Z]+)/) {
            $tagname = $1;
            push (@{$GB{$tagname}}, $line); 
        } else {
            push (@{$GB{$tagname}}, $line); 
        }
    
    }
    
    my ($tag, $version, $GI) = @{$GB{'VERSION'}}[0] =~ m/^([A-Z]+)\s+([^\s]+)\s+([^\s]+)/;
    $self->{'lastload'} = $version;
    $self->{'RECORDS'}->{$version}->{GI} = $GI;
    %{$self->{'RECORDS'}->{$version}->{'RAW'}} = %GB;
    delete ${$self->{'RECORDS'}->{$version}->{'RAW'}}{'VERSION'};
    
}

sub getLastLoaded {
    
    my $self = shift;
    return $self->{'lastload'};
    
}

sub getRecList {
    
    my $self = shift;
    return keys %{$self->{'RECORDS'}};
    
    
}

sub listHeadings {
    
    my $self = shift;
    my $version = shift;
    return keys %{$self->{'RECORDS'}->{$version}};
    
}

sub getRecRef {
    
    my $self = shift;
    my $version = shift;
    return \$self->{'RECORDS'}->{$version};
}

sub processRule {
    
    my $self = shift;
    my $rule = shift;
    my $rulefunct = 'RULE_'.$rule;
    $self->$rulefunct(@_) if exists( &$rulefunct); 
    
}

sub RULE_LOCUS {
    
    my $self = shift;
    my $version = shift;
    $self->{'RECORDS'}->{$version}{'LOCUS'} = genericTextDump(@{$self->{'RECORDS'}->{$version}->{'RAW'}{'LOCUS'}});
    delete($self->{'RECORDS'}->{$version}->{'RAW'}{'LOCUS'});
    
}

sub RULE_ACCESSION {
    
    my $self = shift;
    my $version = shift;
    my ($tag,$value) = split(/\s+/,${$self->{'RECORDS'}->{$version}->{'RAW'}{'ACCESSION'}}[0]);
    
    $self->{'RECORDS'}->{$version}{'ACCESSION'} = $value;
    delete($self->{'RECORDS'}->{$version}->{'RAW'}{'ACCESSION'});
    
}

sub RULE_FEATURES {
    
    my $self = shift;
    my $version = shift;
    my (%features, $tagname, $index, $lasttag);
    
    shift @{$self->{'RECORDS'}->{$version}->{'RAW'}{'FEATURES'}};
    my %block;
    foreach my $line (@{$self->{'RECORDS'}->{$version}->{'RAW'}{'FEATURES'}}) {
        
        if ($line =~ /^\s{5}([A-Za-z]+)/) {
            $tagname = $1;
            my ($lasttag) = keys %block;
            if ( keys %block ) {
                
                my $block = join('',@{$block{$lasttag}});
                my @fields = split(/\s{21}\//, $block);
                my $position = shift @fields;
                $position =~ s/$lasttag|\s//g;
                foreach my $field (@fields) {
                    my ($key,$value);
                    ($key, $value) = $field =~ /^([^=]+)="?([^"]+)"?/;
                    if (! defined $key) {
                        ($key) = $field =~ /^([^=]+)\s?/;
                        if (length $key > 20) {
                            push(@{$self->{'RECORDS'}->{$version}->{$lasttag}{$position}{'note'}}, $value);
                            next;
                        }
                        
                        $value = 'NO DATA';
                        
                        
                    }
                    
                    push(@{$self->{'RECORDS'}->{$version}->{$lasttag}{$position}{$key}}, $value);
                    
                }           
                
                delete $block{$lasttag};
            }
            
            push (@{$block{$tagname}}, $line);
        } else {
            
            push (@{$block{$tagname}}, $line);
        }
    
    }

    delete $self->{'RECORDS'}->{$version}->{'RAW'}{'FEATURES'};
}

sub RULE_ORIGIN {
    
    my $self = shift;
    my $version = shift;
    shift @{$self->{'RECORDS'}->{$version}->{'RAW'}{'ORIGIN'}};
    my $seq = join ("\n",@{$self->{'RECORDS'}->{$version}->{'RAW'}{'ORIGIN'}});
    delete($self->{'RECORDS'}->{$version}->{'RAW'}{'ORIGIN'});
    $seq =~ s/\s|[0-9]|\///g;
    $self->{'RECORDS'}->{$version}{'FULLSEQ'} = $seq;
    
    
}

sub RULE_SOURCE {
    
    my $self = shift;
    my $version = shift;
    shift @{$self->{'RECORDS'}->{$version}->{'RAW'}{'SOURCE'}};
    $self->{'RECORDS'}->{$version}->{'ORGANISM'} = shift @{$self->{'RECORDS'}->{$version}->{'RAW'}{'SOURCE'}};
    $self->{'RECORDS'}->{$version}->{'ORGANISM'} =~ s/\s+ORGANISM\s+//;
    my $taxonomy = join ('',@{$self->{'RECORDS'}->{$version}->{'RAW'}{'SOURCE'}});
    $taxonomy =~ s/\s|\n//g;
    $self->{'RECORDS'}->{$version}->{'TAXONOMY'} = $taxonomy;
    delete $self->{'RECORDS'}->{$version}->{'RAW'}{'SOURCE'};
    
}

sub RULE_DEFINITION {
    my $self = shift;
    my $version = shift;
    my $text = genericTextDump(@{$self->{'RECORDS'}->{$version}->{'RAW'}{'DEFINITION'}});
    $text =~ s/DEFINITION //;
    $self->{'RECORDS'}->{$version}{'DEFINITION'} = $text;
    delete($self->{'RECORDS'}->{$version}->{'RAW'}{'DEFINITION'});
}

sub getOrganism  {
    
    my $self = shift;
    my $version = shift;
    return $self->{'RECORDS'}->{$version}->{'ORGANISM'} if defined $self->{'RECORDS'}->{$version}->{'ORGANISM'};
    $self->processRule('SOURCE',$version);
    return $self->{'RECORDS'}->{$version}->{'ORGANISM'};
    
}

sub getDefinition {
    my $self = shift;
    my $version = shift;
    return $self->{'RECORDS'}->{$version}->{'DEFINITION'} if defined $self->{'RECORDS'}->{$version}->{'DEFINITION'};
    $self->processRule('DEFINITION',$version);
    return $self->{'RECORDS'}->{$version}->{'DEFINITION'};
}

sub getTaxonomy {
    my $self = shift;
    my $version = shift;
    return $self->{'RECORDS'}->{$version}->{'TAXONOMY'} if defined $self->{'RECORDS'}->{$version}->{'TAXONOMY'};
    $self->processRule('SOURCE',$version);
    return $self->{'RECORDS'}->{$version}->{'TAXONOMY'};
}

sub getFullSeq {
    
    my $self = shift;
    my $version = shift;
    return $self->{'RECORDS'}->{$version}{'FULLSEQ'} if defined $self->{'RECORDS'}->{$version}{'FULLSEQ'};
    $self->processRule('ORIGIN',$version);
    return $self->{'RECORDS'}->{$version}{'FULLSEQ'};
    
}

sub getLocus {
    
    my $self = shift;
    my $version = shift;
    return $self->{'RECORDS'}->{$version}{'LOCUS'} if defined $self->{'RECORDS'}->{$version}{'LOCUS'};
    $self->processRule('LOCUS',$version);
    return $self->{'RECORDS'}->{$version}{'LOCUS'};
    
}

sub getAccession {
    
    my $self = shift;
    my $version = shift;
    return $self->{'RECORDS'}->{$version}{'ACCESSION'} if defined $self->{'RECORDS'}->{$version}{'ACCESSION'};
    $self->processRule('ACCESSION',$version);
    return $self->{'RECORDS'}->{$version}{'ACCESSION'};
    
}



sub genericTextDump {
    
    return join("\n",@_);
    
}

sub fillRegionGaps {
    
    my $self = shift;
    my $version = shift;
    my (@start,@stop);    
    
    foreach my $pos ($self->getRegionList($version)) {
        $pos =~ s/<|>//g;
        my ($start, $stop) = $pos =~ /(\d+)\.\.(\d+)/g;
        push (@start,$start);
        push (@stop,$stop);
        
    }
    return if scalar @start != scalar @stop;
    @start = sort {$a <=> $b} @start;
    @stop = sort {$a <=> $b} @stop;
    
    my $end = 0;
    my $regioncount = 1;
    for (my $i = 0; $i < scalar @stop; $i++) {
        
        my $start = $end + 1;
        
        if ($start < $start[$i]) {
            my $regionpos = $start.'..'.($start[$i] - 1);
            my $regionname = 'UNDEFINED-'.$regioncount++;
            $self->{'RECORDS'}->{$version}->{'Region'}{$regionpos}{'region_name'}[0] = $regionname;
        }
        
        $end = $stop[$i];

    }
}

sub hasOrigin {
    
    my $self = shift;
    my $version = shift;
    defined $self->{'RECORDS'}->{$version}->{'RAW'}{'ORIGIN'}? return 1: return 0;
    
}

sub setOrigin {
    
    my $self = shift;
    my $version = shift;
    $self->{'RECORDS'}->{$version}->{'RAW'}{'ORIGIN'} = shift;
    unshift (@{$self->{'RECORDS'}->{$version}->{'RAW'}{'ORIGIN'}},'');

}

sub getRegionList {
    
    my $self = shift;
    my $version = shift;
    
    return keys %{$self->{'RECORDS'}->{$version}->{'Region'}} if !defined $self->{'RECORDS'}->{$version}->{'RAW'}{'FEATURES'};
    $self->processRule('FEATURES',$version);
    return keys %{$self->{'RECORDS'}->{$version}->{'Region'}};
    
}
sub getCDSList {
    
    my $self = shift;
    my $version = shift;
    
    return keys %{$self->{'RECORDS'}->{$version}->{'CDS'}} if !defined $self->{'RECORDS'}->{$version}->{'RAW'}{'FEATURES'};
    $self->processRule('FEATURES',$version);
    return keys %{$self->{'RECORDS'}->{$version}->{'CDS'}};
    
}
sub getCDSData {
    my $self = shift;
    my $version = shift;
    my $location = shift;
    return %{$self->{'RECORDS'}->{$version}->{'CDS'}->{$location}};
    
}
sub getGeneList {
    my $self = shift;
    my $version = shift;
    
    return keys %{$self->{'RECORDS'}->{$version}->{'gene'}} if !defined $self->{'RECORDS'}->{$version}->{'RAW'}{'FEATURES'};
    $self->processRule('FEATURES',$version);
    return keys %{$self->{'RECORDS'}->{$version}->{'gene'}};
}
sub getGeneData {
    my $self = shift;
    my $version = shift;
    my $location = shift;
    return %{$self->{'RECORDS'}->{$version}->{'gene'}->{$location}};
}

sub getmRNAList {
    my $self = shift;
    my $version = shift;
    
    return keys %{$self->{'RECORDS'}->{$version}->{'mRNA'}} if !defined $self->{'RECORDS'}->{$version}->{'RAW'}{'FEATURES'};
    $self->processRule('FEATURES',$version);
    return keys %{$self->{'RECORDS'}->{$version}->{'mRNA'}};
}
sub getmRNAData {
    my $self = shift;
    my $version = shift;
    my $location = shift;
    return %{$self->{'RECORDS'}->{$version}->{'mRNA'}->{$location}};
}

sub getRegionSeq {
    
    my $self = shift;
    my $version = shift;
    my $regionkey = shift;
    return $self->{'RECORDS'}->{$version}->{'Region'}{$regionkey}{'SEQ'} if defined $self->{'RECORDS'}->{$version}->{'Region'}{$regionkey}{'SEQ'};
    $self->parseSeq($version, 'Region', $regionkey);
    return $self->{'RECORDS'}->{$version}->{'Region'}{$regionkey}{'SEQ'};
    
}

sub getGeneSeq {
    my $self = shift;
    my $version = shift;
    my $genekey = shift;
    return $self->{'RECORDS'}->{$version}->{'gene'}{$genekey}{'SEQ'} if defined $self->{'RECORDS'}->{$version}->{'gene'}{$genekey}{'SEQ'};
    $self->parseSeq($version, 'gene', $genekey);
    return $self->{'RECORDS'}->{$version}->{'gene'}{$genekey}{'SEQ'};
}

sub getmRNASeq {
    my $self = shift;
    my $version = shift;
    my $mrnakey = shift;
    return $self->{'RECORDS'}->{$version}->{'mRNA'}{$mrnakey}{'SEQ'} if defined $self->{'RECORDS'}->{$version}->{'mRNA'}{$mrnakey}{'SEQ'};
    $self->parseSeq($version, 'mRNA', $mrnakey);
    return $self->{'RECORDS'}->{$version}->{'mRNA'}{$mrnakey}{'SEQ'};
}

sub getRegionName {
    my $self = shift;
    my $version = shift;
    my $regionkey = shift;
    return $self->{'RECORDS'}->{$version}->{'Region'}{$regionkey}{'region_name'}[0];
}

sub parseSeq {
    
    my $self = shift;
    my $version = shift;
    my $type = shift;
    my $poskey = shift;
    
    my $position = $poskey;
    
    my $completion = 'complete';
    my $complement = $position =~ /complement/;
    my $join = $position =~ /join/;
    my @states = $position =~ /<|>/g;
    $position =~ s/<|>//g;
    
    my %states = ('(complement)' => 'true'
               , '<' => 'Partial on 5` end'
               , '>' => 'Partial on 3` end'
               );
    
    foreach my $code (@states) {
        
        $completion = $states{$code} if $code ne '(complement)';
        
    }
    
    $position =~ s/[^0-9|\.|\,]//g;
    
    
    $self->processRule('ORIGIN',$version) unless defined $self->{'RECORDS'}->{$version}{'FULLSEQ'};
    #my $seq = substr($self->{'RECORDS'}->{$version}{'FULLSEQ'},$start - 1, $end - ($start -1) );
    my $seq = '';
    
    foreach my $segment (split(/,/,$position)) {
        #print $segment,"\n";
        my ($start, $end) = split(/\.\./, $segment );
        $end = $start if ! defined $end;
        $seq .= substr($self->{'RECORDS'}->{$version}{'FULLSEQ'},$start - 1, $end - ($start -1) );
    }    
    
    $seq = _complement_seq($seq) if $complement;
    
    $self->{'RECORDS'}->{$version}->{$type}{$poskey}{'SEQ'} = $seq;
    $self->{'RECORDS'}->{$version}->{$type}{$poskey}{'COMPLETION'} = $completion;
    $self->{'RECORDS'}->{$version}->{$type}{$poskey}{'COMPLEMENT'} = $complement;
    
}

sub _complement_seq {
    
    my $seq = shift;
    $seq = reverse $seq;
    my %comp = ('a' => 't', 't' => 'a', 'g' => 'c', 'c' => 'g','r' => 'y', 'y' => 'r', 's' => 's', 'w' => 'w', 'k' => 'm', 'm' => 'k', 'b' => 'n', 'd' => 'n', 'h' => 'n', 'v' => 'n', 'n' => 'n', '.' => '.', '-' => '-');
    my $comp = '';
    foreach (split(//,$seq)) {
        $comp .= $comp{lc($_)};
    }
    return $comp;
    
=pod
A	Adenine
C	Cytosine
G	Guanine
T (or U)	Thymine (or Uracil)
R	A or G
Y	C or T
S	G or C
W	A or T
K	G or T
M	A or C
B	C or G or T
D	A or G or T
H	A or C or T
V	A or C or G
N	any base
. or -	gap
=cut
    
}


1;
