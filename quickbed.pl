#!/usr/bin/perl -w

# Usage: quickbed.pl [table names...]

# files and directories
my $hg19dir = ".";
my $trackdb = "trackDb";

# default columns to use for BED
my @bedRequired = qw(chrom chromStart chromEnd);
my @bedOptional = qw(name score strand thickStart thickEnd itemRgb blockCount blockSizes blockStarts);
my @bedFieldNames = (@bedRequired, @bedOptional);

# create customized handlers for particular tables (functions that map the table column=>value hashref to individual BED fields)
my %tableHandler =
    ("refGene" => {
	"chromStart" => sub { shift->{"txStart"} },
	"chromEnd" => sub { shift->{"txEnd"} },
	"thickStart" => sub { shift->{"cdsStart"} },
	"thickEnd" => sub { shift->{"cdsEnd"} },
	"blockCount" => sub { shift->{"exonCount"} },
	"blockStarts" => sub { shift->{"exonStarts"} },
	"blockSizes" => sub {
	    my ($rowRef) = @_;
	    my @exonStarts = split (/,/, $rowRef->{"exonStarts"});
	    my @exonEnds = split (/,/, $rowRef->{"exonEnds"});
	    return join (",", map ($exonEnds[$_] - $exonStarts[$_], 0..$#exonEnds));
	},
     },
    );

# Read all .sql files in directory - this was the old way of getting all the tables
# local *DIR;
# opendir DIR, $hg19dir or die "$hg19dir: $!";
# my @sql = map (/^(.*)\.sql$/ ? $1 : (), readdir(DIR));
# closedir DIR;

# Find all the BED files in trackDb.txt.gz
my %trackdbCols = name2column_map ($trackdb);
my %trackdesc;
my @sql;
for_columns ($trackdb,
	     \%trackdbCols,
	     sub {
		 my ($rowRef) = @_;
		 if ($rowRef->{'type'} =~ /^(bed|genePred)\b/) {   # allow bed.* or genePred.* tracks
		     my $tableName = $rowRef->{'tableName'};
		     $trackdesc{$tableName} = $rowRef->{'shortLabel'};
		     push @sql, $tableName;
		 }
	     });

# check that commandline-named BED files are there
for my $table (@ARGV) {
    if (!defined $trackdesc{$table}) {
	warn "Table $table not found in $trackdb\n";
    }
}

# loop over all BED files (or all those named on command line)
local $_;
TABLE: for my $sql (@sql) {   # $sql = name of table (yes, bad choice of variable name)
    # do some checks that we actually want to process this table and all the required files are there
    if (@ARGV && !grep ($_ eq $sql, @ARGV)) {
#	warn "$sql not mentioned in argument list; skipping\n";
	next TABLE;
    }

    my ($sqlfile, $txtgzfile, $bedfile) = map ("$hg19dir/$sql.$_", qw(sql txt.gz bed));

    unless (-e $sqlfile) {
	warn "$sqlfile does not exist; skipping\n";
	next TABLE;
    }

    unless (-e $txtgzfile) {
	warn "$sqlfile does not exist; skipping\n";
	next TABLE;
    }

    if (-e $bedfile) {
	warn "$bedfile already exists; skipping\n";
	next TABLE;
    }

    # create the BED file
    warn "Creating $bedfile from $txtgzfile\n";
    make_bed ($sql, $tableHandler{$sql});
}

# subroutine to crudely parse a .sql table description file and return a map from column names to column indices
sub name2column_map {
    my ($table) = @_;
    my $sqlfile = "$table.sql";

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

    return map (($cols[$_] => $_), 0..$#cols);
}

# subroutine to crudely parse a .txt.gz table dump, create a temporary hashref of column names to values for each row, and apply a given subroutine to this temporary hashref
sub for_columns {
    my ($table, $name2colRef, $func) = @_;

    my $gzip = "gzip -cd $table.txt.gz";
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
	    my %data = map (($_ => $data[$name2colRef->{$_}]), keys %$name2colRef);
	    &$func (\%data);
	    $line = "";
	}
	if (++$lines % 10000 == 0) { warn "(processed $lines/$total lines)\n" }
    }
    close TXT;
}

# wrapper to print a BED file, with an optional hashref mapping BED field names to handler functions that take the table row's colname=>value hashref as input, and return the corresponding BED field as output
sub make_bed {
    my ($table, $handlerRef) = @_;

    # get the column name->index map
    my %name2col = name2column_map ($table);

    # add default handlers
    $handlerRef = {} unless defined $handlerRef;
    for my $bedFieldName (@bedFieldNames) {
	if (!defined($handlerRef->{$bedFieldName})) {
	    if (defined ($name2col{$bedFieldName})) {
		$handlerRef->{$bedFieldName} = sub { shift->{$bedFieldName} };
	    } else {
		$handlerRef->{$bedFieldName} = sub { undef };
	    }
	}
    }

    # create BED file
    my $bedfile = "$table.bed";
    local *BED;
    open BED, ">$bedfile" or die "$bedfile: $!";
    print BED "track name=$table description=\"$trackdesc{$table}\"\n";
    for_columns ($table,
		 \%name2col,
		 sub {
		     my ($tableRowRef) = @_;
		     my @bedFields = map (&{$handlerRef->{$_}} ($tableRowRef), @bedFieldNames);
		     while (@bedFields && !defined($bedFields[$#bedFields])) { pop @bedFields }
		     # could check that required fields are all defined here, I guess...
		     print BED join (" ", map (defined() ? $_ : ".", @bedFields)), "\n";
		 });
    close BED or die "$bedfile: $!";
}
