#! /bin/perl -w

my $filenames=$ARGV[0];
my $output=$ARGV[1];
open OUT2,">$output";
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
	next if $_ =~/^Hugo_Symbol/;
	$gene=(split(/\t/,$_))[0];
	$chr=(split(/\t/,$_))[4];
	$start=(split(/\t/,$_))[5];
	$end=(split(/\t/,$_))[6];
	$ref=(split(/\t/,$_))[11];
	$alt=(split(/\t/,$_))[13];
	$sample=(split(/\t/,$_))[16];
	$info="$gene\t$chr\t$start\t$end\t$ref\t$alt";
	$hash{$info}+=1; #计算info重复的次数
	$data{$info}{$sample}=$_;
}
close I;

foreach my $k (sort keys %data) {
	my $temp;
	foreach my $j (sort keys %{$data{$k}}) {
		 if (exists $hash{$k}>1) {
		 $temp=$data{$k}{$j};
		 #next;
		 }else{
		 	$temp=$data{$k}{$j};#"$k,$hash{$k}","\t",
		 	#next;
		 	#print OUT2 "$data{$k}{$j}\n";
		 }
		 #print OUT2 "$hash{$k}","\t","$data{$k}{$j}\n";
	}
	print OUT2 "$temp\n";

}
close OUT2;

###usage: perl ../Code/redup_HM_isdupsample.pl aml_ohsu_2018_16-01216.maf.txt aml_ohsu_2018_16-01210.maf.txt 








