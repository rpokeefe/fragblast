#!/usr/bin/perl -w
#
# gridgen.pl
# This script takes the ordered output from orderregions.pl and generates heat
# map reports for each query gene for each organism investigated. The output is
# in SVG format and requires an SVG viewer to open. Most web browsers support
# this format now.
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
use SVG;
use NCBI::GenBankEasy;
use File::Slurp;

my @proteins = <global/genbank/*.gp>;
my @organism = qw( physco selag auxeno micro cocco volvox chlamy chlorella ostreo_t ostreo_l);
my %orgs = ('physco' => 'Physcomitrella patens'
            ,'selag' => 'Selaginella moellendorffii'
            ,'auxeno' => 'Auxenochlorella protothecoides'
            ,'cocco' => 'Coccomyxa subellipsoidea'
            ,'micro' => 'Micromonas pusilla'
            ,'volvox' => 'Volvox carteri'
            ,'chlamy' => 'Chlamydomonas reinhardtii'
            ,'chlorella' => 'Chlorella variabilis'
            ,'ostreo_t' => 'Ostreococcus tauri'
            ,'ostreo_l' =>'Ostreococcus lucimarinus');
my @regions;


foreach my $protein (@proteins) {
    undef @regions;
    print "Processing $protein ";
    my $gb = new GenBankEasy;
    my $rec = read_file($protein);
    $gb->loadRec($rec);
    my $version = $gb->getLastLoaded();
    $gb->fillRegionGaps($version);
    
    foreach my $regiontag ($gb->getRegionList($version)) {
    
        my $regionname = $gb->getRegionName($version,$regiontag);
        push @regions,$regiontag.'_'.$regionname;
    
    }
    
    @regions = orderRegions(@regions);
    
    
    my ($outfile) = $protein =~ m/\/([^\/]+)$/;
    my $ordered = $outfile;
    $ordered =~ s/\.gp/\.ordered/;
    $outfile =~ s/\.gp/\.grid\.svg/;
    $outfile = "./grids/".$outfile;
    
    open OUTFILE, ">", $outfile;
    
    my $svg= SVG->new(width=>1600,height=>900);
    
    my $header = $svg->group(
        id    => 'header_'.$version,
    );
    
    $header->text( id=> 'header_'.$version.'_text', x=>800, y=>100, 'text-anchor'=>'middle', 'style'=>{'font-size' => 25}, -cdata=>$version.' '.$gb->getDefinition($version));
    my $midpoint = (75 * scalar @regions) /2;
    my $headeroffset = 800 - $midpoint;
    for (my $i = 0; $i < scalar @regions; $i++) {
        
        $header->rectangle(
            x=>$headeroffset + ($i * 75), y=>150,
            width=>75, height=>50,
            rx=>0, ry=>0,
            id=>'header_'.$regions[$i],
            style=> {
                'stroke' => 'rgb(30,30,30)',
                'stroke-width' => .5,
                'fill' => 'rgb(255,255,255)'
            }
        );
        
        my ($start,$stop,$regname) = $regions[$i] =~ m/(.+)\.\.([^_]+)_(.+)$/;
        $header->text( id=> 'header_'.$regions[$i].'_text', x=>$headeroffset + ($i * 75) + 37, y=>178, 'text-anchor'=>'middle', 'style'=>{'font-size' => 10}, -cdata=>$regname);
        $header->text( id=> 'header_'.$regions[$i].'_start', x=>$headeroffset + ($i * 75) + 37, y=>163, 'text-anchor'=>'middle', 'style'=>{'font-size' => 10}, -cdata=>$start);
        $header->text( id=> 'header_'.$regions[$i].'_stop', x=>$headeroffset + ($i * 75) + 37, y=>193, 'text-anchor'=>'middle', 'style'=>{'font-size' => 10}, -cdata=>$stop);
    }

    my $orgcount = 0;
    foreach my $org (@organism) {
        
        my $orgOrdered = "./$org/ordered/$ordered"; 
        renderFile($svg, $orgOrdered, $org,$orgcount);
        print ".";
        $orgcount++;
    }
    
    my $datestring = localtime();
    my $footer = $svg->group(
        id    => 'footer_'.$version,
    );
    
    $footer->text( id=> 'footer_'.$version.'_text', x=>800, y=>700, 'text-anchor'=>'middle', 'style'=>{'font-size' => 12}, -cdata=>"FragBlast Report   $datestring");
    $footer->text( id=> 'footer_'.$version.'_name', x=>800, y=>715, 'text-anchor'=>'middle', 'style'=>{'font-size' => 12}, -cdata=>"Ronald O'Keefe");
    
    print OUTFILE $svg->xmlify;
    close OUTFILE;
    print "\n";
}


sub renderFile {
    my $svg = shift;
    my $filename = shift;
    my $orgname = shift;
    my $orgnum = shift;
    
    
       
    my @lines = read_file($filename);
    my $header = shift @lines;
    chomp $header;
    my @headerRegions = split(/,/,$header);
    shift @headerRegions;
    
    my @map;
    my $ordercount = 0;
    
    for ( my $i = 0; $i < scalar @regions; $i++) {
        $regions[$i] =~ s/<|>//g;
        if ($regions[$i] eq $headerRegions[$ordercount]) {
            
            $map[$i] = $ordercount++;
            
        } else {
            $map[$i] = undef;
        }
        
    }
    
    
    
    my @grid;
    
    for (my $i = 0; $i < 10; $i++ ) {
        last if ! defined $lines[$i];
        my $line = $lines[$i];
        my @frags = split(/,/,$line);
        my ($version, $length) = splice @frags,0,2;
        
        my $limit = scalar @map;
        
        for ( my $j = 0; $j < $limit; $j++)  {
            
            next unless defined $map[$j];
            
            my $eval = getBestEval($frags[$map[$j]]);
            $grid[$i][$j] = $eval; 
        
        }
        
    }
    
    my $gridwidth = scalar @map * 15;
    
    my $gridy = (int($orgnum / 5) * 200) + 300;
    my $gridx = ($orgnum % 5  * 300) + 50 + 150 - ($gridwidth/2);
    
    drawThumbPrint($svg, $orgname,\@grid, $gridx, $gridy);

}

sub drawThumbPrint {
    my $svg = shift;
    my $organism = shift;
    my @grid = @{shift @_};
    my $posx = shift;
    my $posy = shift;
    
    my $gridcol = 0;
    if (defined $grid[0]) { $gridcol = scalar @{$grid[0]}}
    my $gridrow = scalar @grid;
    #my @startpos = (100,100);
    
    my $y=$svg->group(
        id    => $organism,
    );
    
    
    my $textx = $posx + $gridcol * 7.5;
    my $texty = $posy - 10;
    $y->text( id=> $organism.'_text', x=>$textx, y=>$texty, 'text-anchor'=>'middle', 'style'=>{'font-size' => 15}, -cdata=>$orgs{$organism});
    return if $gridcol == 0;
    for (my $row = 0; $row < $gridrow ; $row++) {
        
        my $starty = $posy + $row*15;
        
        for (my $col = 0; $col < $gridcol; $col++) {
            
            my $eval = $grid[$row][$col];
            my $startx = $posx + $col*15;
            my $id = $organism.'_'.$startx.'_'.$starty;
            drawBox($y, $startx, $starty, $id, $eval);
            
        }
        
    }

}

sub getBestEval {
    
    my $frag = shift;
    my @evals = sort {$a <=> $b } $frag =~ m/:([^:]+):/g;
    return shift @evals;
    
}

sub drawBox {
    
    my $y = shift;
    my $startx = shift;
    my $starty = shift;
    my $id = shift;
    my $eval = shift;
    
    $y->rectangle(
        x=>$startx, y=>$starty,
        width=>15, height=>15,
        rx=>0, ry=>0,
        id=>$id,
        style=> {
            'stroke' => 'rgb(30,30,30)',
            'stroke-width' => .5,
            'fill' => getColor($eval)
        }
    );    
    
}

sub getColor {
    
    my $eval = shift;
    return 'rgb(0,0,0)' unless $eval;
    $eval = $eval**-1;
    my $max = 1e-100**-1;
    
    my @colorgrid = ([0,undef,255]
                     ,[0,255,undef]
                     ,[undef,255,0]
                     ,[255,undef,0]
                     ,[255,0,0]);
    
    my $rating = log10($eval)/log10($max) * 1024;
    $rating = 1024 if $rating > 1024;
    $rating = 0 if $rating < 0;
    
    my $offset = $rating % 256;
    my $row = int ($rating/256);
    
    for (my $i = 0; $i < 3; $i++) {
        
        if (! defined ${$colorgrid[$row]}[$i] ) {
            
            $offset = 256 - $offset if $row % 2;
            ${$colorgrid[$row]}[$i] = $offset;
    
        }
        
    }
    my $rgb = join(',',@{$colorgrid[$row]});
    return "rgb($rgb)";
    
}

sub log10 {
    my $n = shift;
    return log($n)/log(10);
}


sub orderRegions {
    
    my %regions;
    #shift @_;
    foreach (@_) {
        
        unless (/FULL/) {
            
            chomp;
            $_ =~ s/<|>//g;
            my ($start,$finish,$name) = /(\d+)\.\.(\d+)_(.+)/;
            die "$start already exists: $name" if exists $regions{$start};
            $regions{$start} = $_;
            
        }
        
    }
    
    return map {$regions{$_}} sort{ $a <=> $b } keys %regions;
    
}


