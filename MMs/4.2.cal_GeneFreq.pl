#! /bin/perl -w

my $path=$ARGV[0];
my (@filename,$f);
## method1:  calculate all sample from each tumor, which cost a lot time
opendir DIR,"$path";
@filenames=grep{/Multi_SameGene.txt$/}readdir DIR;  
closedir(DIR);

open OUT,">Pan.Multi_SameGene.Freq.txt";
print OUT "#Cancer\tGene\tPts\tMultiPts\tMutPts\tTotalPts\tFreqinMulti\tFreqinMut\tFreqinTotal\n";

my $Mutpath=$ARGV[1];  ##具体突变的位置
#step1:read data(file)
my $gene;
for (my $i=0; $i<@filenames; $i++){
	open I, "$path\/$filenames[$i]"; 
	my $cancer=(split(/\./,$filenames[$i]))[0];
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
		###构建患者的列表
		 push @Temp,$patient;
		 #print @Pts,"\n";
	}
	close I;
   @Pts=grep {++$hash{$_}==1} @Temp; #代码块{ ++$count{$_} == 1 }作为无重复元素的判断条件，并通过grep得到一个满足无重复元素的子数组;去除冗余元素	 ,具有基因多突变模式的患者 
   @Uniq=grep {++$hash{$_}==1} @Uniqtemp;#去除冗余元素	 ,所有患者
#print @Uniqtemp,"\n";
#####基因频率
	my %count;
	foreach  (@Genelist) {$count{ $_ }++;}

	foreach my $a (keys %count) {
		my $freq=$count{$a}/($#Pts+1);  #### 具有某基因多突变模式的患者\/具有基因多突变模式的所有患者
		my $freqinTotal=$count{$a}/($#Uniq+1+$#Pts+1);  #### 具有某基因多突变模式的患者\/所有患者
        ##提取具有该突变的患者
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
		for (my $i=0; $i<@maffile; $i++){
			open MUT, "$folder\/$maffile[$i]";
			readline MUT;
			my $index=0;
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
			}
			next if ($index==0);
			push @PtsMut,$maffile[$i];
		 }
		close MUT;
		#print $cancer,"\t",$a,"\t",$count{$a},"\t",$#PtsMut+1,"\t",$folder,"\n";
		my $freqinMUT=$count{$a}/($#PtsMut+1); ###具有某基因多突变模式的患者\/具有该基因的突变的所有患者
		print OUT $cancer,"\t",$a,"\t",$count{$a},"\t",$#Pts+1,"\t",$#PtsMut+1,"\t",$#Uniq+1+$#Pts+1,"\t",$freq,"\t",$freqinMUT,"\t",$freqinTotal,"\n";

		}

}

###usage method1: perl ../../Code/cal_GeneFreq.pl /p200/liuxin_group/moshl/Temp/TCGA/Results/Multi/temp  /p200/liuxin_group/moshl/Temp/TCGA        

