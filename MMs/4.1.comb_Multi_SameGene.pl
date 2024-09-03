#! /bin/perl -w

my $path=$ARGV[0];
my (@filename,$f);
opendir DIR,"$path";
@filenames=grep{/Multi_SameGene.txt$/}readdir DIR;
closedir(DIR);

open OUT1,">Pan.Multi_SameGene.txt";
print OUT1 "#Cancer\tSample\tGene\tCount\tType\tIndel\tINS\tInfo\n";

#step1:read data(file)
my $gene;
for (my $i=0; $i<@filenames; $i++)
{
	open I, "$path\/$filenames[$i]"; 
	my $cancer=(split(/\./,$filenames[$i]))[0];
	readline I;
	my @Genelist=();
	my %hash;
	my (@Pts,$patient,@Temp,@Uniq,@Uniqtemp);
	while (<I>)
	{
		chomp;
		$gene=(split(/\t/,$_))[1];
		$patient=(split(/\t/,$_))[0];
		push @Uniqtemp,$patient;
		next if !$gene;
		print OUT1 $cancer,"\t",$_,"\n";
	}
	close I;
}

###usage perl ../../Code/comb_Multi_SameGene.pl /p200/liuxin_group/moshl/Temp/TCGA/Results/Multi/temp   
