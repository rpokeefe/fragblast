#!/usr/bin/perl -w
#
# bubblegen.pl
# This script generates the SVG bubble charts using the ordered files.
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
use SVG;

my @organism = qw( micro volvox auxeno cocco ostreo_t ostreo_l chlamy chlorella selag physco );

foreach my $organism (@organism) {

    my @ordered = <$organism/ordered/*.ordered>;
    foreach my $filename (@ordered) {
        print "processing [$organism] $filename\n";
        my $infile = my $outfile = $filename;
        $outfile =~ s/ordered/bubble/g;
        $outfile .= ".svg";
        drawBlast($infile,$outfile);
        
    }
}

sub drawBlast {
    my $file = shift;
    my $outfile = shift;
    my $basewidth = 1200;
    my $ypos = 75;
    
    my @rec = read_file($file);
    my ($geneversion) = $file =~ /ordered\/(.+)___Arab/;
    
    my $docheight = 150 + (10 * 100);
    
    my @gene = split(/,/,shift @rec);
    shift @gene;
    my $last = $gene[$#gene];
    my ($genelength) = $last =~ /\.\.(\d+)_/;
    return unless $genelength;
    
    my $max = $genelength;
    my $reccounter = 0;
    foreach my $rec (@rec) {
        
        last if $reccounter++ == 10;
        my @list = split(/,/,$rec);
        $max = $list[1] if $list[1] > $max;
    
    }
    
    my $scale = $basewidth/$max;
    
    my $svg = SVG->new(width=>$basewidth,height=>$docheight);
    
    # draw query gene
    my $y = $svg->group(
        id    => 'group_y',
        style => {
            'font-family' => 'arial'
        }
    );
    
    addLabel($y,$geneversion,$ypos);
    addBaseLine($y,$geneversion,$genelength * $scale,$ypos);
    
    my %usednames;
    my @domains;
    
    my $limit = scalar @gene;
    
    for  (my $i = 0; $i < $limit; $i++) {
        
        my $frag = $gene[$i];
        
        $frag =~ s/(UND)EFINED\-(\d+)/$1$2/;
        my ($start,$end,$tag) = $frag =~ /^(\d+)\.\.(\d+)_(.+)$/;
        push(@domains,$tag);
        #next if exists $nested{$i};
        
        
        if (exists $usednames{$tag}) {
            $usednames{$tag}++;
        } else {
            $usednames{$tag} = 1;
        }
        my $tagid = $tag."|".$usednames{$tag};
        
        my $boxfill = 'rgb(200,255,200)';
        my $boxstroke = 'rgb(0,255,0)';
        my $fontcolor = 'black';
        my $opacity = 0.8;
        
        if ($tag =~ /UND/) {
            $boxfill = 'black';
            $boxstroke = 'black';
            $fontcolor = 'white';
            $opacity = 0.5;
        }
        
        
        $start = 0 if $start == 1;
        my $startadj = $start * $scale;
        my $endadj = $end * $scale;
        
        addBubble($y,$tagid,$tag,$start,$end,$endadj -$startadj,$startadj, $ypos, $boxfill, $boxstroke, $fontcolor);
        
    }
    
    my $counter = 0;
    foreach my $hit (@rec) {
        
        last if $counter++ == 10;
        chomp $hit;
        my @hits = split(/,/,$hit);
        my $version = shift @hits;
        my $genelength = shift @hits;
        my $limit = scalar @hits;
        
        $ypos += 100;
        
        my %usedid;
        
        my %fill = ( DELTA => 'rgb(255,200,200)', PSSM => 'rgb(200,255,200)', BLASTP => 'rgb(255,255,200)', PSIBLAST => 'rgb(200,200,255)');
        my %stroke = ( DELTA => 'rgb(255,0,0)', PSSM => 'rgb(0,255,0)', BLASTP => 'rgb(255,102,0)', PSIBLAST => 'rgb(0,0,255)');
        my %font = ( DELTA => 'black', PSSM => 'black', BLASTP => 'black', PSIBLAST => 'black');
        
        my $group = $svg->group(
            id    => 'group_'.$version,
            style => {
                'font-family' => 'arial'
            }
        );
        
        addLabel($group,$version,$ypos);
        addBaseLine($group,$version,$genelength*$scale,$ypos);
        
        for (my $i = 0; $i < $limit; $i++ ) {
            
            my @types = split(/\|/, $hits[$i]);
            
            foreach (@types) {
                
                my ($type,$escore,$start,$stop) = split(/:|\.\./, $_);
                my $id = $version.'_'.$type.'_'.$domains[$i];
                
                if (exists $usedid{$id}) {
                    
                    $usedid{$id}++;
                    
                } else {
                    $usedid{$id} = 1;
                }
                
                $id .= "|".$usedid{$id};
                
                addBubble($group, $id ,$domains[$i], $start, $stop, ($stop - $start) * $scale, $start * $scale, $ypos, $fill{$type}, $stroke{$type}, $font{$type});
                
                
            }
            
        }
    
    }
    
    
    open OUTPUT, ">", $outfile or die $!;
    print OUTPUT $svg->xmlify;
    close OUTPUT;
    
}

sub addLabel {
    my $group = shift;
    my $label = shift;
    my $ypos = shift;
    
    $group->text( id=>$label,
             x=>10, y=>$ypos - 35,
             #'text-anchor'=>'middle',
             style => {
                'font-size' => 16,
                'fill' => 'black'
             },
             -cdata=>$label
             );
}

sub addBaseLine {
    my $group = shift; #group reference
    my $version = shift;
    my $width = shift;
    my $ypos = shift;
    
    $group->line(
        id=>$version.'_base',
        x1=>0, y1=>$ypos,
        x2=>$width, y2=>$ypos,
        style => {
            'stroke' => 'rgb(0,0,0)',
            'stroke-width' => 5,
            'opacity' => 0.3
        }
    );
    
}

sub addBubble {
    
    my $group = shift;
    my $bubbleid = shift;
    my $bubblename = shift;
    my $start = shift;
    my $stop = shift;
    my $length = shift;
    my $xpos = shift;
    my $ypos = shift;
    my $fill = shift;
    my $stroke = shift;
    my $font = shift;
    my $opacity = 0.8;
    
    $group->rectangle(
        x=>$xpos, y=>$ypos - 25,
        width=>$length, height=>50,
        rx=>10, ry=>10,
        id=>$bubbleid,
        style=> {
            'stroke' => $stroke,
            'stroke-width' => 1,
            'fill' => $fill,
            'fill-opacity' => $opacity
        }
    );
    
    my $textx = $xpos + ($length/2);
    $group->text( id=>$bubbleid.'_text',
             x=>$textx, y=>$ypos + 2, #78
             'text-anchor'=>'middle',
             style => {
                'font-size' => 4,
                'fill' => $font
             },
             -cdata=>$bubblename
             );
    
    $group->text( id=>$bubbleid.'_start',
             x=>$textx, y=>$ypos - 18, #58,
             'text-anchor'=>'middle',
             style => {
                'font-size' => 4,
                'fill' => $font
             },
             -cdata=>$start
             );
    $group->text( id=>$bubbleid.'_end',
             x=>$textx, y=>$ypos + 22, #98,
             'text-anchor'=>'middle',
             style => {
                'font-size' => 4,
                'fill' => $font
             },
             -cdata=>$stop
             );

}
