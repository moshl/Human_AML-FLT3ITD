#! /bin/perl -w
use POSIX;
#1.ѡ��������������: ԭ��������
my $samplefile=$ARGV[0];
my (%hash,$idx);
my (%Sample,%Patient,$pts);
open  Samp,"$samplefile";
while (<Samp>) {
	next if $_ !~/^TCGA/;
	chomp $_;
	$idx=(split /\t/,$_)[1];
	$pts=(split /\t/,$_)[0];
    next if ($idx !~/01$/ && $idx !~/03$/ && $idx !~/09$/ ); # && $idx !~/05$/  #05�������ж�Ӧ��01������Ϊ�˼���ļ���ԣ���������
	$hash{$idx}=1;
	$Patient{$pts}+=1;
	$Sample{$pts}{$idx}=1;
}

#2.ѡ�����ͻ�����ݵ�����
my $mutlist=$ARGV[1];
open I, $mutlist; 
readline I;
my @samplename;
my %Patients;
while (<I>)
{
	next if $_ !~/^case_list_ids/;
	chomp $_;
	my @tmp=split '\t', (split /\: /,$_)[1];
	foreach my $id (@tmp) {
		next if (not exists $hash{$id});
		push @samplename,$id;
	}
	
}
#3.�ҳ��ظ�����
foreach my $k (sort keys %Sample) {
	my $info="";
	foreach my $j (sort keys %{$Sample{$k}}) {
		if (grep {$_ eq $j} @samplename){
		if ($Patient{$k}>1) {
			$info.=$j;
		}else{
				$info=$j;
			}
		}
	}
	#next if (length($info)<2);
	print $k,"\t",$Patient{$k},"\t",$info,"\n"; ##�������Ҫ��
}

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


##run dictionary: eg. /p200/liuxin_group/moshl/Temp/TCGA/Data/LAML
##usage perl ../../Code/SplitMAF_TCGA_bySample.pl ./data_bcr_clinical_data_sample.txt ./case_lists/cases_sequenced.txt ./data_mutations_mskcc.txt > IsDupSampleInfo.txt
##cd ../UCEC/ && perl ../../Code/SplitMAF_TCGA_bySample.pl data_bcr_clinical_data_sample.txt case_lists/cases_sequenced.txt data_mutations_mskcc.txt 33 34 0 > IsDupSampleInfo.txt

