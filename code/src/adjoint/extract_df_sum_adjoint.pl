#======================================================================================================================!
#
# extract_df_sum_adjoint.pl
#
#======================================================================================================================!
#
# finish_to_gen_adjoint.pl unconditionally deletes every generated file whose
# name matches /mpi_(back|diff)/ (tap/m_mpi_back.f90, tap/m_mpi_diff.f90) --
# historically because Tapenade's own (wrong) differentiation of com_var_r/
# com_dof lives there and must be discarded in favour of the hand-written
# com_var_r_back/com_dof_back (m_mpi_back.f90 in src/adjoint).
#
# Since df_sum_r/df_sum_i were renamed away from mpi_sum_r/mpi_sum_i so that
# Tapenade would stop misrecognizing them as opaque MPI primitives (see the
# comment on df_sum_r in src/common/m_mpi.f90), Tapenade now ALSO correctly
# auto-generates DF_SUM_R_BACK/DF_SUM_R_DIFF (calling MPI_ALLREDUCE_FWD/BWD/D
# from the ADFirstAidKit) inside that same doomed m_mpi_back.f90/m_mpi_diff.f90
# -- and would be deleted right along with the unwanted com_var_r content.
#
# This script runs BEFORE finish_to_gen_adjoint.pl and rescues just the
# DF_SUM_R_BACK / DF_SUM_R_DIFF subroutines into their own file, named so it
# does not match the /mpi_(back|diff)/ deletion pattern and survives. Being
# pulled out of their enclosing MODULE M_MPI_BACK/M_MPI_DIFF, they need their
# own USE/INCLUDE/declarations added back for the module-level variables
# (val_tmp_r, code) and Tapenade-introduced temporaries (val_tmp_r_back/diff)
# they reference.
#
# (DF_SUM_I has no _BACK/_DIFF: mpi_sum_i is only ever used on plain integer
# counters that are inactive w.r.t. the control vector, so Tapenade
# correctly determines it needs no derivative.)

use strict;
use warnings;

sub extract_subroutine {
    my ($text, $name) = @_;
    if ( $text =~ /(^\s*SUBROUTINE\s+\Q$name\E\b.*?^\s*END\s+SUBROUTINE\s+\Q$name\E.*?\n)/ims ) {
        return $1;
    }
    return undef;
}

my $back_file = 'm_mpi_back.f90';
my $diff_file = 'm_mpi_diff.f90';

if ( -f $back_file ) {

    open( my $in, '<', $back_file ) or die "cannot open $back_file: $!";
    local $/;
    my $content = <$in>;
    close $in;

    my $body = extract_subroutine( $content, 'DF_SUM_R_BACK' );

    if ( $body ) {
        $body =~ s/(\bIMPLICIT NONE\b)/USE m_common\n    USE m_mpi, only: val_tmp_r, code\n    $1\n    INCLUDE 'admpif.h'\n    REAL(rp) :: val_tmp_r_back/;
        open( my $out, '>', 'df_sum_reduce_back.f90' ) or die $!;
        print $out "! Rescued from $back_file by extract_df_sum_adjoint.pl (see comment there):\n";
        print $out "! finish_to_gen_adjoint.pl deletes any *mpi_back*/*mpi_diff* file wholesale,\n";
        print $out "! which would otherwise discard this Tapenade-generated, correct adjoint too.\n";
        print $out "$body\n";
        close $out;
        print "extract_df_sum_adjoint.pl: rescued DF_SUM_R_BACK -> df_sum_reduce_back.f90\n";
    } else {
        print "extract_df_sum_adjoint.pl: WARNING DF_SUM_R_BACK not found in $back_file -- Tapenade's output shape may have changed, adjust this script.\n";
    }
}

if ( -f $diff_file ) {

    open( my $in, '<', $diff_file ) or die "cannot open $diff_file: $!";
    local $/;
    my $content = <$in>;
    close $in;

    my $body = extract_subroutine( $content, 'DF_SUM_R_DIFF' );

    if ( $body ) {
        # DF_WAIT_ALL is a module procedure of m_mpi (mangled __m_mpi_MOD_df_wait_all),
        # not an external one (df_wait_all_) -- must USE it explicitly now that this
        # subroutine no longer lives inside MODULE M_MPI_DIFF (which itself USEd m_mpi).
        $body =~ s/(\bIMPLICIT NONE\b)/USE m_common\n    USE m_mpi, only: val_tmp_r, code, df_wait_all\n    $1\n    INCLUDE 'admpif.h'\n    REAL(rp) :: val_tmp_r_diff/;
        open( my $out, '>', 'df_sum_reduce_diff.f90' ) or die $!;
        print $out "! Rescued from $diff_file by extract_df_sum_adjoint.pl -- see df_sum_reduce_back.f90\n";
        print $out "$body\n";
        close $out;
        print "extract_df_sum_adjoint.pl: rescued DF_SUM_R_DIFF -> df_sum_reduce_diff.f90\n";
    } else {
        print "extract_df_sum_adjoint.pl: WARNING DF_SUM_R_DIFF not found in $diff_file -- Tapenade's output shape may have changed, adjust this script.\n";
    }
}
