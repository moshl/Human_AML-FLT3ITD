#! /bin/perl -w
my $samefile=$ARGV[0];
my (%hash,$id,$gene);
open  Same,"$samefile";
readline Same;
while (<Same>) {
	$id=(split /\./,((split /\t/,$_)[0]))[0];
	#print $id,"\n";
	$gene=(split /\t/,$_)[1];
	#print $gene,"\n";
	next if !$gene;
    $hash{$id}{$gene}=1;
}
close Same;
#print $samplename[1],"\n";

my $mutfile=$ARGV[1];
#print $header,"\n";
my $header=`head -1 $mutfile `;
my $name=$ARGV[2];
open OUT,">$name.Multi_SameGene.MUTInfo.txt";
print OUT $header;
open I, $mutfile; 
readline I;
my @Class=('Missense_Mutation','Nonsense_Mutation','Nonstop_Mutation','Splice_Site','Splice_Region','Frame_Shift_Del','Frame_Shift_Ins','In_Frame_Del','In_Frame_Ins');
my ($cls,$num_cls);
while (<I>)
{
	chomp $_;
	$cls=(split(/\t/,$_))[9];
	$num_cls = grep /$cls/i, @Class; #print $num_cls,"\n";
	next if ($num_cls <1);  #next if ($_ =~/Silent/ | );
	my $hugo=(split(/\t/,$_))[0];
	my $sample=(split(/\t/,$_))[16];
	#print $sample,"\n";
	if (exists $hash{$sample}{$hugo}){
		print OUT $_,"\n";
	}	
	
}
close I;

##usage perl ../../Code/selectSample_MMIS.pl ./temp/LAML.Multi_SameGene.txt ../../Data/LAML/data_mutations_mskcc.txt  LAML







