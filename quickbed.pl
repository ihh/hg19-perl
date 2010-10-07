#!/usr/bin/perl -w

# Usage: quickbed.pl [table names...]

my $hg19dir = ".";

my @bedRequired = qw(chrom chromStart chromEnd);
my @bedOptional = qw(name score strand);

my $trackdb = "trackDb";

# Read all .sql files in directory - this was the old way of getting all the tables
# local *DIR;
# opendir DIR, $hg19dir or die "$hg19dir: $!";
# my @sql = map (/^(.*)\.sql$/ ? $1 : (), readdir(DIR));
# closedir DIR;

# Find all the BED files in trackDb.txt.gz
my @trackdbCols = name2column_list ("$trackdb.sql", qw(tableName shortLabel type));
my %trackdesc;
my @sql;
for_columns ("$trackdb.txt.gz",
	     \@trackdbCols,
	     sub {
		 my ($tableName, $shortLabel, $type) = @_;
		 if ($type =~ /^bed\b/) {
		     $trackdesc{$tableName} = $shortLabel;
		     push @sql, $tableName;
		 }
	     });

# check that commandline-named BED files are there
for my $file (@ARGV) {
    if (!defined $trackdesc{$file}) {
	warn "Table $file not found in $trackdb\n";
    }
}

# loop over all BED files (or all those named on command line)
local $_;
for my $sql (@sql) {
    my ($sqlfile, $txtgzfile, $bedfile) = map ("$hg19dir/$sql.$_", qw(sql txt.gz bed));

    if (@ARGV && !grep ($_ eq $sql, @ARGV)) {
#	warn "$sql not mentioned in argument list; skipping\n";
	next;
    }

    unless (-e $txtgzfile) {
	warn "$txtgzfile does not exist; skipping\n";
	next;
    }

    if (-e $bedfile) {
	warn "$bedfile already exists; skipping\n";
	next;
    }

    my %name2col = name2column_map ($sqlfile, @bedRequired, @bedOptional);
    my @missing = grep(!defined($name2col{$_}), @bedRequired);
    if (@missing) {
	warn "$sql does not have required column(s): @missing\n";
	next;
    }

    my @order = map ($name2col{$_}, @bedRequired, @bedOptional);
    while (!defined($order[$#order])) { pop @order }

    warn "Converting $txtgzfile to $bedfile\n";

    open BED, ">$bedfile" or die "$bedfile: $!";
    print BED "track name=$sql description=\"$trackdesc{$sql}\"\n";
    for_columns ($txtgzfile,
		 \@order,
		 sub { print BED join (" ", map (defined() ? $_ : ".", @_)), "\n" });
    close BED or die "$bedfile: $!";
}

# subroutine to crudely parse a .sql table description file and return a map from column names to column indices
sub name2column_map {
    my ($sqlfile, @colNames) = @_;

    my @cols;
    local *SQL;
    local $_;
    open SQL, "<$sqlfile" or die "$sqlfile: $!";
    while (<SQL>) { last if /CREATE TABLE/ }
    while (<SQL>) {
	last if /^\)/;
	if (/^\s*\`(\S+)\`/) { push @cols, $1 }
    }
    close SQL;

    return map ($cols[$_] =~ /^(@{[join("|",@colNames)]})$/
		? ($1 => $_)
		: (),
		0..$#cols);
}

# wrapper subroutine to return the values of name2column_map as an ordered list
sub name2column_list {
    my ($sqlfile, @colNames) = @_;
    my %n2c = name2column_map ($sqlfile, @colNames);
    return map ($n2c{$_}, @colNames);
}

# subroutine to crudely parse a .txt.gz table dump, select the indexed columns, and apply a given subroutine to each row so selected
sub for_columns {
    my ($txtgzfile, $orderRef, $func) = @_;

    my $gzip = "gzip -cd $txtgzfile";
    my $total = `$gzip | wc -l` + 0;

    local *TXT;
    local $_;
    open TXT, "$gzip|" or die "$gzip: $!";
    my $lines = 0;
    my $line = "";
    while (<TXT>) {
	chomp;
	if (/\\$/) {
	    chop;
	    $line .= $_;
	} else {
	    $line .= $_;
	    my @data = split /\t/, $line;
	    &$func (map (defined() && defined($data[$_])
			 ? $data[$_]
			 : undef,
			 @$orderRef));
	    $line = "";
	}
	if (++$lines % 10000 == 0) { warn "(processed $lines/$total lines)\n" }
    }
    close TXT;
}

