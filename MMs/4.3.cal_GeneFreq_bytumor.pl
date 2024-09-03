#! /bin/perl -w

## method2: calculation sample from one type of tumor,
my $path=$ARGV[0];
my $filenames=$ARGV[1];
my $Mutpath=$ARGV[2];  ##����ͻ���λ��
#step1:read data(file)
my $gene;
open I, "$path/$filenames"; 
my $cancer=(split(/\./,$filenames))[0];
print $cancer,"\n";
open OUT,">$cancer.Multi_SameGene.Freq.txt";
print OUT "#Cancer\tGene\tInDel\tSNV\tInDelpts\tSNVpts\tPts\tMultiPts\tMutPts\tTotalPts\tFreqinMulti\tFreqinMut\tFreqinTotal\n";

readline I;
my @Genelist=();
my %hash;
my (@Pts,$patient,@Temp,@Uniq,@Uniqtemp);
while (<I>){
	chomp;
	$gene=(split(/\t/,$_))[1];
	$patient=(split(/\t/,$_))[0];
	push @Uniqtemp,$patient;
	next if !$gene;
	#print "$gene\n";
	push @Genelist,$gene;
	#$patient=(split(/\t/,$_))[0];
	###�������ߵ��б�
	 push @Temp,$patient;
	 #print @Pts,"\n";
}
close I;
@Pts=grep {++$hash{$_}==1} @Temp; #�����{ ++$count{$_} == 1 }��Ϊ���ظ�Ԫ�ص��ж���������ͨ��grep�õ�һ���������ظ�Ԫ�ص�������;ȥ������Ԫ��	 ,���л����ͻ��ģʽ�Ļ��� 
@Uniq=grep {++$hash{$_}==1} @Uniqtemp;#ȥ������Ԫ��	 ,���л���
#print @Uniqtemp,"\n";
#####����Ƶ��
my %count;
foreach  (@Genelist) {$count{ $_ }++;}

foreach my $a (keys %count) {
	my $freq=$count{$a}/($#Pts+1);  #### ����ĳ�����ͻ��ģʽ�Ļ���\/���л����ͻ��ģʽ�����л���
	my $freqinTotal=$count{$a}/($#Uniq+1+$#Pts+1);  #### ����ĳ�����ͻ��ģʽ�Ļ���\/���л���
    ##��ȡ���и�ͻ��Ļ���
	my (@PtsMut,@MUTAll);
	my $folder;
	if ($cancer=~/^HM/) {
		my $tmp=(split(/HM/,$cancer))[1];
		$folder="$Mutpath\/HM/$tmp/"; 
	}else{
		$folder="$Mutpath\/Data/$cancer/"; 
	}
	opendir DIR,"$folder";
	my @maffile=grep{/maf.uniq.txt$/}readdir DIR;
	my @Class=('Missense_Mutation','Nonsense_Mutation','Nonstop_Mutation','Splice_Site','Splice_Region','Frame_Shift_Del','Frame_Shift_Ins','In_Frame_Del','In_Frame_Ins');
	my @Variant_type=('SNP','INS','DEL');
	my $SNV=0;
	my $InDel=0;
	my $SNVpts=0;
	my $InDelpts=0;
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
	        next if ($chr !~ /^[0-9]+/); #�ų�X��Y��MT�ϵ�ͻ��
	        next if ($num_cls<1 || $num_var<1);  #next if ($_ =~/Silent/ | );				
			my $hugo=(split(/\t/,$_))[0];
			next if ($a !~ $hugo);
			$index++;
			my $type=(split(/\t/,$_))[10];
			if ($type=~/NP$/) {
				$SNV++;
				$snvidx++
			}else{
				$InDel++;
				$indelidx++;
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
	my $freqinMUT=$count{$a}/($#PtsMut+1); ###����ĳ�����ͻ��ģʽ�Ļ���\/���иû����ͻ������л���
	print OUT $cancer,"\t",$a,"\t",$InDel,"\t",$SNV,"\t",$InDelpts,"\t",$SNVpts,"\t",$count{$a},"\t",$#Pts+1,"\t",$#PtsMut+1,"\t",$#Uniq+1+$#Pts+1,"\t",$freq,"\t",$freqinMUT,"\t",$freqinTotal,"\n";

}


###usage method2: perl ../../Code/cal_GeneFreq_bytumor.pl /p200/liuxin_group/moshl/Temp/TCGA/Results/Multi/temp LAML.Multi_SameGene.txt  /p200/liuxin_group/moshl/Temp/TCGA       

