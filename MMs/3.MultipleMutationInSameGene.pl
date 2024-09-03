#! /bin/perl -w

my $path=$ARGV[0];
my (@filename,$f);
opendir DIR,"$path";
@filenames=grep{/maf.uniq.txt$/}readdir DIR;
#foreach $f(@filenames){
#print "$path\/$f \n";
##system "grep cd $f |wc -l \n";
#}
closedir(DIR);

my $hname=$ARGV[1];
open OUT,">$hname.Multi_SameGene.txt";
print OUT "#Sample\tGene\tCount\tType\tIndel\tINS\tInfo\n";

open OUT1,">$hname.Freq.txt";
print OUT1 "#Cancer\tMulitNum\tTotal\tFreq\n";


open OUT2,">$hname.Multi_SameGene.MUTInfo.txt";
my $header=`head -1 $path/$filenames[1] `;
print OUT2 $header;

#step1:read data(file)
my %hash=();
my ($gene,$cls,$num_cls,$var,$num_var,$chr);
my $num=0;
my @Class=('Missense_Mutation','Nonsense_Mutation','Nonstop_Mutation','Splice_Site','Splice_Region','Frame_Shift_Del','Frame_Shift_Ins','In_Frame_Del','In_Frame_Ins');
my @Variant_type=('SNP','INS','DEL');
##for循环，0ne-by-one
for (my $i=0; $i<@filenames; $i++){
	open I, "$path\/$filenames[$i]"; #print $filenames[$i],"\n";
	readline I;
	my @Genelist=();
	my @multi_array=();
	my @Mut;
	while (<I>)
	{
		chomp;
		next if (!$_); #跳过空行
		##1.select key mutations
		$cls=(split(/\t/,$_))[9]; 
		$var = (split(/\t/,$_))[10];
		$chr = (split(/\t/,$_))[4]; #print $cls,"\t",$var,"\t",$chr,"\n";
		$num_cls = grep /$cls/i, @Class; #print $num_cls,"\n";
		$num_var = grep /$var/i, @Variant_type; #print $num_var,"\n";
		next if ($chr !~ /^[0-9]+/); #排除X、Y和MT上的突变
		next if ($num_cls<1 || $num_var<1);  #next if ($_ =~/Silent/ | );
		$gene=(split(/\t/,$_))[0];
		#print $gene,"\n";
		next if ($gene=~/^\./);
		push @Genelist,$gene;
		$hash{$gene}{$filenames[$i]}=$_;
		push @Mut,$_;
	}
	close I;

	my %Temp;  #p.s该数组不可与后面计算元素出现次数用的数组一致
    my @dup=grep { ++$Temp{ $_ } > 1; } @Genelist;  ##这是判断语句，给出来的结果是最新添加进的元素是否重复，是，输出基因，否，输出空
    #print @dup,"\n";

	my %count;
	foreach  (@Genelist) {
		$count{ $_ }++;  ##每个基因出现的次数
	}
	if (@dup) {
		$num++;
		my $target;
		foreach  (keys %count) {
		   next if ($count{$_}<=1);
		   print OUT $filenames[$i],"\t", $_,"\t",$count{$_};
		   $target=$_;
		   my @temp;
		   my $info;
		   foreach $item (@Mut) {
			   my $g=(split(/\t/,$item))[0];
			   next if ($g !~ $target || length($g) > length($target) || length($g) < length($target));
			   print OUT2 $item,"\n";
			   my $type=(split(/\t/,$item))[10];
			   my $HGVSp_Short=(split(/\t/,$item))[39];
			   my $ref_count=(split(/\t/,$item))[33];
			   my $alt_count=(split(/\t/,$item))[34];
			   $info.="\t$type\t$HGVSp_Short\t$ref_count\t$alt_count";
			   #print OUT "\t" ,$type,"\t",$ref_count,"\t",$alt_count;
			   push @temp,$type;
		   }
		   my $InDel=0;
		   my $INS=0;
		   foreach $a (@temp) {
			   next if($a =~/SNP$/); #只考虑SNP，排除DNP，TNP
			   $InDel++;
			   if($a =~/INS/){
				   $INS++;
			   }
		   }
		   if($InDel==0){
			   print OUT "\t",1,"\t",$InDel,"\t",$INS,$info,"\n"; ##only SNV
		   }elsif ($InDel==$count{$_}) {
			   print OUT "\t",2,"\t",$InDel,"\t",$INS,$info,"\n";  ##only InDel
		   }else{
			   print OUT "\t",3,"\t",$InDel,"\t",$INS,$info,"\n";  ## mixed, SNV/InDel
		   }
		  
		}
	}else{
		$num=$num;
		print OUT  $filenames[$i],"\n";
	}
    
	## 计算 Type and Diff

}


my $total=@filenames;
my $freq=$num/$total;
print OUT1 "$hname\t$num\t$total\t$freq\n";

###usage: perl ../../../Code/MultipleMutationInSameGene.pl /p200/liuxin_group/moshl/Temp/TCGA/Data/LAML LAML








