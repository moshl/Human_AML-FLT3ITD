#! /bin/perl -w

my $filenames=$ARGV[0];
open OUT2,">data_mutations_mskcc_rmdup.txt";
my $header=`head -1 $filenames `;
print OUT2 $header;

#step1:read data(file)
my %hash=();
my ($gene,$chr,$start,$end,$ref,$alt,$sample);
my %data;

open I, "$filenames"; 
readline I;
while (<I>)
{
	chomp;
	##1.select key mutations
	$gene=(split(/\t/,$_))[0];
	$chr=(split(/\t/,$_))[4];
	$start=(split(/\t/,$_))[5];
	$end=(split(/\t/,$_))[6];
	$ref=(split(/\t/,$_))[11];
	$alt=(split(/\t/,$_))[12];
	$sample=(split(/\t/,$_))[16];
	$info="$gene\t$chr\t$start\t$end\t$ref\t$alt\t$sample";
	$hash{$info}+=1; #计算info重复的次数
	$caller=(split(/\t/,$_))[32];
	$data{$info}{$caller}=$_;
}
close I;

foreach my $k (sort keys %hash) {
	if (exists $data{$k}{mutect}) {
		print OUT2 "$data{$k}{mutect}\n"; #"$k,$hash{$k}","\t",
		next;
		}elsif (exists $data{$k}{pindel}) {
			print OUT2 "$data{$k}{pindel}\n"; #"$k,$hash{$k}","\t",
			next;
		}else{
			print OUT2 "$data{$k}{varscan}\n"; #"$k,$hash{$k}","\t",
			next;
		}

}
close OUT2;

###usage: perl ../Code/redup_HM.pl data_mutations_mskcc.txt








