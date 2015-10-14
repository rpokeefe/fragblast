package BlastLocalv2;
#############################################################################
# BlastLocalv2.pm
# This package provides a PERL interface to the Blast command line tool provided
# by NCBI.
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

sub newProfile {
    
    my $self = shift;
    my $profile = shift;
    $self->{'PROFILE'}->{$profile}->{'PROGRAM'} = 'psiblast';
    $self->{'PROFILE'}->{$profile}->{'PARAMETERS'}->{'-db'} = 'chlorella';
    $self->{'PROFILE'}->{$profile}->{'PARAMETERS'}->{'-matrix'} = 'BLOSUM45';
    $self->{'PROFILE'}->{$profile}->{'PARAMETERS'}->{'-gapopen'} = 15;
    $self->{'PROFILE'}->{$profile}->{'PARAMETERS'}->{'-gapextend'} = 2;
    $self->{'PROFILE'}->{$profile}->{'PARAMETERS'}->{'-num_alignments'} = 1000;
    $self->{'PROFILE'}->{$profile}->{'PARAMETERS'}->{'-num_descriptions'} = 1000;
    
}

sub deleteProfile {
    my $self = shift;
    my $profile = shift;
    
    delete $self->{'PROFILE'}->{$profile} if exists $self->{'PROFILE'}->{$profile}; 
    
}
sub set_db {
    
    my $self = shift;
    my $profile = shift;
    my $db = shift;
    $self->{'PROFILE'}->{$profile}->{'PARAMETERS'}->{'-db'} = $db;
    
}

sub getProfiles {
    
    my $self = shift;
    return keys %{$self->{'PROFILE'}};
    
}

sub setProgram {
    my $self = shift;
    my $profile = shift;
    my $program = shift;
    if ($program =~ /blastn|blastx|blastp|tblastn|tblastx|psiblast|rpsblast|rpsblastn|deltablast/) {
        $self->{'PROFILE'}->{$profile}->{'PROGRAM'} = $program;
    } else {
        die "Program: $program does not have valid option";
    }
}

sub set_out_pssm {
    my $self = shift;
    my $profile = shift;
    my $file = shift;
    
    if ( $self->{'PROFILE'}->{$profile}->{'PROGRAM'} eq 'psiblast' ) {
        $self->{'PROFILE'}->{$profile}->{'PARAMETERS'}->{'-out_pssm'} = $file;
        $self->{'PROFILE'}->{$profile}->{'PARAMETERS'}->{'-out_ascii_pssm'} = $file.'.ascii';
    }
    
}

sub set_comp_based_stats {
    
    my $self = shift;
    my $profile = shift;
    
    $self->{'PROFILE'}->{$profile}->{'PARAMETERS'}->{'-comp_based_stats'} = shift;
    
}

sub unset_out_pssm {
    
    my $self = shift;
    my $profile = shift;
    
    if ( defined $self->{'PROFILE'}->{$profile}->{'PARAMETERS'}->{'-out_pssm'} ) {
        delete $self->{'PROFILE'}->{$profile}->{'PARAMETERS'}->{'-out_pssm'};
        delete $self->{'PROFILE'}->{$profile}->{'PARAMETERS'}->{'-out_ascii_pssm'};
    }
    
    
}

sub set_in_pssm {
    
    my $self = shift;
    my $profile = shift;
    my $file = shift;
    
    if ( $self->{'PROFILE'}->{$profile}->{'PROGRAM'} eq 'psiblast' ) {
        $self->{'PROFILE'}->{$profile}->{'PARAMETERS'}->{'-in_pssm'} = $file; 
    }
    
}

sub set_num_iterations {
    
    my $self = shift;
    my $profile = shift;
    my $iterations = shift;
    
    if ( $self->{'PROFILE'}->{$profile}->{'PROGRAM'} eq 'psiblast' ) {
        $self->{'PROFILE'}->{$profile}->{'PARAMETERS'}->{'-num_iterations'} = $iterations;
    }
    
}

sub useBLOSUM62 {
    
    my $self = shift;
    my $profile = shift;
    $self->{'PROFILE'}->{$profile}->{'PARAMETERS'}->{'-matrix'} = 'BLOSUM62';
    $self->{'PROFILE'}->{$profile}->{'PARAMETERS'}->{'-gapopen'} = 11;
    $self->{'PROFILE'}->{$profile}->{'PARAMETERS'}->{'-gapextend'} = 1;
    
}

sub unset_matrix {
    
    my $self = shift;
    my $profile = shift;
    
    if ( defined $self->{'PROFILE'}->{$profile}->{'PARAMETERS'}->{'-matrix'} ) {
    
        delete $self->{'PROFILE'}->{$profile}->{'PARAMETERS'}->{'-matrix'};
        delete $self->{'PROFILE'}->{$profile}->{'PARAMETERS'}->{'-gapopen'};
        delete $self->{'PROFILE'}->{$profile}->{'PARAMETERS'}->{'-gapextend'};
    
    }
    
}

sub set_query {
    
    my $self = shift;
    my $profile = shift;
    my $seq = shift;
    $self->{'PROFILE'}->{$profile}->{'PARAMETERS'}->{'-query'} = $seq;
    
}

sub unset_query {
    
    my $self = shift;
    my $profile = shift;
    
    delete $self->{'PROFILE'}->{$profile}->{'PARAMETERS'}->{'-query'} if defined $self->{'PROFILE'}->{$profile}->{'PARAMETERS'}->{'-query'};
    
}

sub set_out {
    
    my $self = shift;
    my $profile = shift;
    $self->{'PROFILE'}->{$profile}->{'PARAMETERS'}->{'-out'} = shift;
    
}

sub buildCommand {
    
    my $self = shift;
    my $profile = shift;
    my %params = %{$self->{'PROFILE'}->{$profile}->{'PARAMETERS'}};
    $self->{'PROFILE'}->{$profile}->{'COMMAND'} = join (' ',$self->{'PROFILE'}->{$profile}->{'PROGRAM'} ,map { $_.' '.$params{$_}} keys %params);
    return $self->{'PROFILE'}->{$profile}->{'COMMAND'};
    
}

sub execute {
    
    my $self = shift;
    my $profile = shift;
    
    local $| = 1;
    return system($self->buildCommand($profile));
    
}

1;
