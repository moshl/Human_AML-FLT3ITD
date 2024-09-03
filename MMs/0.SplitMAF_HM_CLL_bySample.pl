#! /bin/perl -w
my (%hash,%Patient,%Sample);
my $idx;
my $tempfile=$ARGV[0];
open IN,"$tempfile";
while (<IN>) {
	next if $_ !~/treatment naive/;
	#print $_,"\n";
	chomp $_;
	$idx=(split /\t/,$_)[0];
	$hash{$idx}=1;
}
close (IN);

my $samplefile=$ARGV[1];
my (@samplename,$id);
my $pts;
open  Samp,"$samplefile";
while (<Samp>) {
	chomp $_;
	next if $_ !~/^CLL/;
	$id=(split /\t/,$_)[1];
	$pts=(split /\t/,$_)[0];
	if (exists $hash{$pts}){
	#print $pts,"\n";
	$Patient{$pts}+=1;
	$Sample{$pts}{$id}=1;
	push @samplename,$id;
	}
}
close (Samp);

foreach my $k (sort keys %Sample) {
	my $info;
	foreach my $j (sort keys %{$Sample{$k}}) {
		if ($Patient{$k}>1) {
			$info.=$j;
		}else{
				$info=$j;
			}
	}
	print $k,"\t",$Patient{$k},"\t",$info,"\n"; ##�������Ҫ��
}
#print $samplename[1],"\n";

#my %hash=();
#my $mutfile=$ARGV[2];
##print $header,"\n";
#my $header=`head -1 $mutfile `;
##print  $header,"\n"; #OUT
#
#for (my $i=0; $i<@samplename; $i++){
#	open OUT,">$samplename[$i].maf.txt"; 
#	print OUT $header;
#	open I, $mutfile; 
#	readline I;
#	while (<I>)
#	{
#		chomp $_;
#		next if $_ !~ $samplename[$i];
#		print OUT $_,"\n";
#		
#	}
#	close I;
#}

##3.�ָ��ļ���ɾ���ظ�ͻ��
my $mutfile=$ARGV[2];
my $header=`head -1 $mutfile `;
my $ref_col=$ARGV[3];
my $alt_col=$ARGV[4];
my $cls=$ARGV[5];
for (my $i=0; $i<@samplename; $i++){
	open OUT,">$samplename[$i].maf.txt"; 
	print OUT $header;
	open I, $mutfile; 
	readline I;
	my (%mut,$type);
	while (<I>)
	{
		chomp $_;
		next if $_ !~ $samplename[$i];
		my @snv=split(/\t/,$_);
		#ɾ��ͬһλ�õĶ��ͻ��
		my $chr=$snv[4];
		my $start=$snv[5];	   
		my $end=$snv[6];
		my $ref=$snv[11];
		my $alt=$snv[13];
		#$type="$chr:$start:$end:$ref:$alt";
		$type="$chr:$start"; #1.�ų�ͬһ��λ���ж����ͬ�����͵�ͻ��, ����ͬ����ͻ��ѡ�������ߵ�����ͻ���¼
		#print $type,"\n";
		my $alt_count=$snv[$alt_col];
		#��ֵ��Ϊ�˺���depth�ļ���		
		if(!$alt_count || !( $alt_count =~ /^[0-9]+/) ){
			$alt_count=0;
		}
		my $ref_count=$snv[$ref_col];
		#print $ref_count, "\n";
		if(!$ref_count || !( $ref_count =~ /^[0-9]+/) ){
			$ref_count=0;
		}
		
		if($cls==1){
			$ref_count=$ref_count-$alt_count;   ##total reads - alt reads
		}elsif ($cls==2){
			my $tmp=floor $ref_count*$alt_count/100; ##total reads * VAF
			$alt_count=$tmp;
			$ref_count=$ref_count-$alt_count;
		}
		##�滻��Ϊ�����εļ���
        $snv[33]=$ref_count;
		$snv[34]=$alt_count;				
		#print $ref_count,"\t",$alt_count,"\n";
		my $depth=$ref_count+$alt_count; 
		#print $ref_count,"\t",$depth,"\n";
		print OUT join("\t",@snv),"\n";
		
		if(!exists $mut{$type}){
		 $mut{$type}=join("\t",@snv);	#	$mut{$type}=$_	
		}else{
			my $RD=(split(/\t/,$mut{$type}))[33];
		    if(!$RD){
		    	$RD=0;
		    }
		    my $AD=(split(/\t/,$mut{$type}))[34];
		    if(!$AD){
		    	$AD=0;
		    }
		    my $DP=$RD+$AD; 
			#print $RD,"\t",$DP,"\n";
			if($DP<$depth){
				delete $mut{$type};
				$mut{$type}=join("\t",@snv);	#	$mut{$type}=$_
			}
		}
		
	}
	close I;
	
	my (%rmdup,$gene,@rest);
	foreach my $k (sort keys %mut){
		 my $gene=(split(/\t/,$mut{$k}))[0];
		 my $s1=(split(/\t/,$mut{$k}))[5];
		 my $r1=(split(/\t/,$mut{$k}))[33];
		 if(!$r1 || !( $r1 =~ /^[0-9]+/) ){
			$r1=0;
		 }
		 my $a1=(split(/\t/,$mut{$k}))[34];
		 if(!$a1 || !( $a1 =~ /^[0-9]+/) ){
			$a1=0;
		 }
		 my $d1=$r1+$a1;
		 if(!exists $rmdup{$gene}){
		 $rmdup{$gene}=$mut{$k};			
		 }else{
			my $first=(split(/\t/,$rmdup{$gene}))[5];						
			if(abs($s1-$first)<50){ #2.�ų��ٽ�λ�õ�ͻ��
				my $RD1=(split(/\t/,$rmdup{$gene}))[33];
		        if(!$RD1){
		        	$RD1=0;
		        }
		        my $AD1=(split(/\t/,$rmdup{$gene}))[34];
		        if(!$AD1){
		        	$AD1=0;
		        }
		        #$alt_count=undef($alt_count);
		        my $DP1=$RD1+$AD1; 
		 	   #print $RD,"\t",$DP,"\n";
		 	   if($DP1<$d1){
		 	   	delete $rmdup{$gene};
		 	   	$rmdup{$gene}=$mut{$k};	
		 	   }
			}else{
				#print $rmdup{$gene},"\n";
				#print $mut{$k},"\n";
				#print $s1,"\t",$first,"\t",abs($s1-$first), "\n";
				push @rest,$rmdup{$gene},"\n";
				$rmdup{$gene}=$mut{$k};
				#print $rmdup{$gene},"\n";
			}
		 }	  
	}
	
    open OUT2,">$samplename[$i].maf.uniq.txt"; 
	print OUT2 $header;
	foreach my $g (sort keys %rmdup){
		 print OUT2 $rmdup{$g},"\n";
	}
	print OUT2 @rest,"\n";
	
	
}


##run dictionary: eg. /p200/liuxin_group/moshl/Temp/TCGA/HM/CLL
##usage perl ../../Code/SplitMAF_HM_CLL_bySample.pl ./data_clinical_patient.txt ./data_clinical_sample.txt ./data_mutations_mskcc.txt 
