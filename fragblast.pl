#!/usr/bin/perl -w
#
# fragblast.pl
# This script parses the arabidopsis protein query genbank files,
# selects the protein regions and blast each regaion against the
# listed organism blast databases. Four different blast methods
# used for each region.
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
use NCBI::GenBankEasy;
use NCBI::BlastLocalv2;
use File::Slurp;

# A list of abbreviated organism names. These names are used as blast database names
# as well as folder names.

my @organism = qw( micro volvox auxeno cocco ostreo_t ostreo_l chlamy chlorella selag physco );

# Get a list of Genbank files that will be used for queries.
my @files = <global/genbank/*.gp>;

# process each subject organism.
foreach my $organism (@organism) {

    # process each query file.
    foreach my $file (@files) {
    
        print "[$organism] Processing:  ",$file,"\n";
        
        # create genbank file object and create regions.
        my $gb = new GenBankEasy;
        my $rec = read_file($file);
        $gb->loadRec($rec);
        my $version = $gb->getLastLoaded();
        $gb->fillRegionGaps($version);
        my ($id) = $file =~ m/genbank\/(.+)\.gp/;
        
        #perform blasts on the full query.
        my $blast = new BlastLocalv2;    
        my $fullprofile = $id.'___FULL';
        my $seq = $gb->getFullSeq($version);
        blastp($blast,$fullprofile,$seq,$organism);
        psiblast($blast,$fullprofile,$seq,$organism);
        pssmblast($blast,$fullprofile,$seq,$organism);
        deltablast($blast,$fullprofile,$seq,$organism);
        
        # perform blasts on all regions of query protein. 
        foreach my $regiontag ($gb->getRegionList($version)) {
        
            my $regionname = $gb->getRegionName($version,$regiontag);
            my $profile = $id.'___'.$regiontag.'_'.$regionname;
            $profile =~ s/<|>//g;
            my $seq = $gb->getRegionSeq($version,$regiontag);
            #add all blast blocks here
            blastp($blast,$profile,$seq,$organism);
            psiblast($blast,$profile,$seq,$organism);
            pssmblast($blast,$profile,$seq,$organism);
            deltablast($blast,$profile,$seq,$organism);
        
        }
    }
}

# subroutine executes deltablast
sub deltablast { 
    
    my $blast = shift;
    my $profile = shift;
    my $seq = shift;
    my $organism = shift;
    
    print "[$organism] executing deltablast on $profile\n";
    
    # sets input and output files for blast.
    my $query = './global/fasta/'.$profile.'.fasta';
    my $blastout = "./$organism/blastout/".$profile.'____DELTA.txt';
    
    # if input fasta file hasn't been created, create it.
    if (! -e $query) {
        open FASTA, '>', $query or die $!;
        print FASTA '>', $profile,"\n",$seq;
        close FASTA;
    }
    
    # delete any existing profiles with the same identifier and create a new one.
    $blast->deleteProfile($profile);
    $blast->newProfile($profile);
    # set parameters for blast.
    $blast->setProgram($profile,'deltablast');
    $blast->unset_matrix($profile);
    $blast->set_query($profile,$query);
    $blast->set_db($profile,$organism);
    $blast->set_num_iterations($profile,3);
    # set outfile and execute blast.
    $blast->set_out($profile,$blastout);
    $blast->execute($profile);
    
}

# subrouting executes protein blasts. Structurally similiar to sub deltablast.
sub blastp {
    
    my $blast = shift;
    my $profile = shift;
    my $seq = shift;
    my $organism = shift;
    
    print "[$organism] executing blastp on $profile\n";
    
    my $query = './global/fasta/'.$profile.'.fasta';
    my $blastout = "./$organism/blastout/".$profile.'____BLASTP.txt';
    
    if (! -e $query) {
        open FASTA, '>', $query or die $!;
        print FASTA '>', $profile,"\n",$seq;
        close FASTA;
    }
    
    #build search chlorella
    $blast->deleteProfile($profile);
    $blast->newProfile($profile);
    $blast->setProgram($profile,'blastp');
    $blast->set_query($profile,$query);
    $blast->set_db($profile,$organism);
    $blast->set_out($profile,$blastout);
    $blast->execute($profile);
    
}

# subrouting executes PSI blasts. Structurally similiar to sub deltablast.
sub psiblast {
    
    my $blast = shift;
    my $profile = shift;
    my $seq = shift;
    my $organism = shift;
    
    print "[$organism] executing phiblast $profile\n";
    
    my $query = './global/fasta/'.$profile.'.fasta';
    my $blastout = "./$organism/blastout/".$profile.'____PSIBLAST.txt';
    
    if (! -e $query) {
        open FASTA, '>', $query or die $!;
        print FASTA '>', $profile,"\n",$seq;
        close FASTA;
    }
    
    #build search chlorella
    $blast->deleteProfile($profile);
    $blast->newProfile($profile);
    $blast->setProgram($profile,'psiblast');
    $blast->set_comp_based_stats($profile,1);
    $blast->set_query($profile,$query);
    $blast->set_db($profile,$organism);
    $blast->set_num_iterations($profile,3);
    $blast->set_out($profile,$blastout);
    $blast->execute($profile);
    
}

#subroutine executes a psi blast with a generated pssm file.
sub pssmblast { 
    
    my $blast = shift;
    my $profile = shift;
    my $seq = shift;
    my $organism = shift;
    
    print "[$organism] executing pssmblast $profile\n";
    
    my $query = './global/fasta/'.$profile.'.fasta';
    my $blastout = "./$organism/blastout/".$profile.'____PSSM.txt';
    my $pssmblastout = './global/pssmblastout/'.$profile.'.txt';
    my $pssm = './global/pssm/'.$profile.'.pssm';
    
    if (! -e $query) {
        open FASTA, '>', $query or die $!;
        print FASTA '>', $profile,"\n",$seq;
        close FASTA;
    }
    
    $blast->deleteProfile($profile);
    $blast->newProfile($profile);
    $blast->setProgram($profile,'psiblast');
    $blast->set_comp_based_stats($profile,1);
    
    # build a pssm for the fragment if it does not exist.
    if ( ! -e $pssm) {
        print "$profile -BUILDING PSSM\n";
        # weakphyte is a blastdb created to handle the pssm creation.
        # it contains organism databases that had a weaker homology to
        # phytochrome. (auxeno, cocco, chlamy, ostreo_l, ostreo_t, volvox)
        $blast->set_db($profile,'weakphyte');
        $blast->set_query($profile,$query);
        $blast->set_out($profile,$pssmblastout);
        $blast->set_num_iterations($profile,2);
        $blast->set_out_pssm($profile,$pssm);
        $blast->execute($profile);
    }
    
    # set pssm file for query and blast with pssm.
    print "$profile -SCANNING TARGET\n";
    if (-e $pssm) { #sometimes a pssm isn't generated
        $blast->unset_query($profile);
        $blast->set_in_pssm($profile,$pssm);
    } else {
        print "$profile -NO PSSM FOUND: using query file\n";
    }
    $blast->set_db($profile,$organism);
    $blast->unset_out_pssm($profile);
    $blast->set_num_iterations($profile,1);
    $blast->set_out($profile,$blastout);
    $blast->execute($profile);
}
