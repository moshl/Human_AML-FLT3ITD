#! /bin/perl -w

## background
my $cancer=$ARGV[0];
my $Mutpath=$ARGV[1];  ##具体突变的位置

print $cancer,"\n";
open OUT,">$cancer.TotalGene.Freq2.txt";
#print OUT "#Cancer\tGene\tInDel\tSNV\tInDelpts\tSNVpts\tMutPts\n";

#step1:read data(file)
#构建突变基因列表
my @Class=('Missense_Mutation','Nonsense_Mutation','Nonstop_Mutation','Splice_Site','Splice_Region','Frame_Shift_Del','Frame_Shift_Ins','In_Frame_Del','In_Frame_Ins');
my @Variant_type=('SNP','INS','DEL');

my $folder;
if ($cancer=~/^HM/) {
	my $tmp=(split(/HM/,$cancer))[1];
	$folder="$Mutpath\/HM/$tmp/"; 
}else{
	$folder="$Mutpath\/Data/$cancer/"; 
}
opendir DIR,"$folder";


#已经完成的基因
my $ok=$ARGV[2];
open I, $ok; 
readline I;
my (%hash,$g);
while (<I>)
{
	chomp $_;
	next if (!$_); #跳过空行
	$g=(split(/\t/,$_))[1];
	$hash{$g}=1;
    #print $g,"\n";	
}

my @mutlist=grep{/mutations_mskcc.txt$/}readdir DIR;
#print @mutlist;
#####识别所有突变基因
my $gene;
open I, "$folder\/$mutlist[0]"; 
readline I;
my %count;
while (<I>)
{
	my $a=grep /$_/i, @Class; 
	next if ($a <1);
	chomp $_;
	$gene=(split(/\t/,$_))[0];
	next if (exists $hash{$gene});
	#print $gene,"\n";
	$count{$gene}=1;
	
	
}
#####基因频率
foreach my $a (keys %count) {
    ##提取具有该突变的患者
	my (@PtsMut,@MUTAll);
	opendir DIR,"$folder";
	my @maffile=grep{/maf.uniq.txt$/}readdir DIR;
	my $SNV=0;
	my $InDel=0;
	my $SNVpts=0;
	my $InDelpts=0;
	my $len_info=",";
	for (my $i=0; $i<@maffile; $i++){
		open MUT, "$folder\/$maffile[$i]";
		readline MUT;
		my $index=0;
		my $snvidx=0;
		my $indelidx=0;
		while (<MUT>) {
			chomp;
			next if (!$_);
			##1.select key mutations
			my $cls=(split(/\t/,$_))[9];
			my $num_cls = grep /$cls/i, @Class; #print $num_cls,"\n";
	        my $var = (split(/\t/,$_))[10];
	        my $num_var = grep /$var/i, @Variant_type; #print $num_var,"\n";
	        my $chr = (split(/\t/,$_))[4]; #print $cls,"\t",$var,"\t",$chr,"\n";
	        next if ($chr !~ /^[0-9]+/); #排除X、Y和MT上的突变
	        next if ($num_cls<1 || $num_var<1);  #next if ($_ =~/Silent/ | );				
			my $hugo=(split(/\t/,$_))[0];
			next if ($a !~ $hugo);
			$index++;
			my $type=(split(/\t/,$_))[10];
			if ($type=~/SNP$/) {
				$SNV++;
				$snvidx++
			}else{
				$InDel++;
				$indelidx++;
				##2. Indel长度
			    my $ref=(split(/\t/,$_))[11];
			    my $alt=(split(/\t/,$_))[13];
				my $len;
				if ($type=~/INS$/) {
					$len=length($alt);
				}else{
					$len=0-length($ref);
				}
				$len_info.=$len."\,";
			}
			

			
		}
		close MUT;
		next if ($index==0);
		push @PtsMut,$maffile[$i];
		if ($snvidx>0) {
			$SNVpts++;
		}
		if ($indelidx>0) {
			$InDelpts++;
		}
	}
	
	#print $cancer,"\t",$a,"\t",$count{$a},"\t",$#PtsMut+1,"\t",$folder,"\n";
	print OUT $cancer,"\t",$a,"\t",$InDel,"\t",$SNV,"\t",$InDelpts,"\t",$SNVpts,"\t",$#PtsMut+1,"\t",$len_info,"\n";

}


###usage: perl ../../Code/cal_GeneFreq_bytumor_total.pl LAML /p200/liuxin_group/moshl/Temp/TCGA
