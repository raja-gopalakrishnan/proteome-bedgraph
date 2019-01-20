# C is minus
# W is plus
library(stringr)
options(scipen=100)

#Input the fasta file that has sequences of all proteins
fasta=read.delim2(paste(getwd(),"/orf_trans.fasta",sep=""), header=FALSE)
fasta=apply(fasta,2,as.character)

warning_geneID=NULL
amino_lines=apply(fasta,1,function(x){!identical(as.character(substring(x,1,1)),">")})
header_lines=which(!amino_lines)

#Select proteins encoded by the nuclear genome - exclude mitochondrial ones
nucleus_headers=which(!(as.character(substring(fasta[header_lines,1],2,2))!="Y"))

bedgraph=matrix(nrow=0,ncol=6)

for(i in nucleus_headers)
{
  #Extract details for the protein from the header line
  chr_coord=regexpr("Chr",fasta[header_lines[i],1])[1]+4
  strand=substring(fasta[header_lines[i],1],8,8)
  
  start=regexpr("from",fasta[header_lines[i],1])[1]+5
  end=regexpr(", G",fasta[header_lines[i],1])[1]-1
  
  chromosome=paste("chr",substring(fasta[header_lines[i],1],chr_coord,(start-7)),sep="")
  coord=substring(fasta[header_lines[i],1],start,end)
  first_space=regexpr(" ",fasta[header_lines[i],1])[1]
  geneID=substring(fasta[header_lines[i],1],2,first_space-1)
  
  exon_coord=as.data.frame(strsplit(coord,","))
  exon_coord=t(as.data.frame(apply(exon_coord,1,function(x){as.data.frame(strsplit(x,"-"))})))
  
  #Extract amino acid sequence and convert to bedgraph for coding region on crick strand
  if(strand=="C")
  {
    if(nrow(exon_coord)!=1)
    {
      exon_coord=as.data.frame(apply(exon_coord,2,rev),stringsAsFactors=FALSE)
    }
    
    #Extract amino acid sequene for the protein
    amino_acid=character()
    aa_rows=(fasta[(header_lines[i]+1):(header_lines[i+1]-1),])
    for(q in 1:length(aa_rows))
    {
      amino_acid=str_c(amino_acid,aa_rows[q])
    }
    
    col2=NULL
    col3=NULL
    aa=NULL
    diff2=0
    
    for(k in 1:nrow(exon_coord))
    {
      diff=((as.numeric(exon_coord[(k),1])-as.numeric(exon_coord[(k),2])+1)%%3)
      diff1=3-diff2
      diff2=(diff-diff1)%%3
      diff3=diff2+(1.5*(diff2-1)*(diff2-2))
      
      quotient=(as.numeric(exon_coord[(k),1])-as.numeric(exon_coord[(k),2]))>2
      if(quotient)
      {
        new_col2=c(exon_coord[k,1],seq((as.numeric(exon_coord[k,1])-diff1),(as.numeric(exon_coord[k,2])+diff3-1),-3))
        new_col3=c(seq((as.numeric(exon_coord[k,1])-diff1+1),(as.numeric(exon_coord[k,2])+diff2),-3),exon_coord[k,2])
        new_col2=new_col2[!duplicated(new_col2)]
        new_col3=new_col3[!duplicated(new_col3)]
        col2=c(col2,new_col2)
        col3=c(col3,new_col3)
      }else
      {
        new_col2=exon_coord[k,1]
        new_col3=exon_coord[k,2]
        col2=c(col2,new_col2)
        col3=c(col3,new_col3)
      }
      if(diff1==3)
      {
        if(k<3) temp=substring(amino_acid,(length(aa)+1),(length(aa)+length(new_col2)))
        else temp=substring(amino_acid,(length(aa)+3-k),(length(aa)+length(new_col2)+2-k))
        aa=c(aa,as.character(strsplit(temp,"")[[1]]))
      }
      else
      {
        if(k<3) temp=substring(amino_acid,(length(aa)),(length(aa)+length(new_col2)-1))
        else temp=substring(amino_acid,(length(aa)+2-k),(length(aa)+length(new_col2)+1-k))
        aa=c(aa,as.character(strsplit(temp,"")[[1]]))
      }
    }
    bedgraph=rbind(bedgraph,cbind(chromosome,col3,col2,"-",aa,geneID))
  }
  
  #Extract amino acid sequence and convert to bedgraph for coding region on watson strand
  if(strand=="W")
  {
    amino_acid=character()
    aa_rows=(fasta[(header_lines[i]+1):(header_lines[i+1]-1),])
    for(q in 1:length(aa_rows))
    {
      amino_acid=str_c(amino_acid,aa_rows[q])
    }
    
    col2=NULL
    col3=NULL
    aa=NULL
    diff2=0
    
    for(k in 1:nrow(exon_coord))
    {
      diff=((as.numeric(exon_coord[(k),2])-as.numeric(exon_coord[(k),1])+1)%%3)
      diff1=3-diff2
      diff2=(diff-diff1)%%3
      diff3=diff2+(1.5*(diff2-1)*(diff2-2))
      
      quotient=(as.numeric(exon_coord[(k),2])-as.numeric(exon_coord[(k),1]))>2
      if(quotient)
      {
        new_col2=c(exon_coord[k,1],seq((as.numeric(exon_coord[k,1])+diff1),(as.numeric(exon_coord[k,2])-diff3+1),3))
        new_col3=c(seq((as.numeric(exon_coord[k,1])+diff1-1),(as.numeric(exon_coord[k,2])-diff2),3),exon_coord[k,2])
        new_col2=new_col2[!duplicated(new_col2)]
        new_col3=new_col3[!duplicated(new_col3)]
        col2=c(col2,new_col2)
        col3=c(col3,new_col3)
      }else
      {
        new_col2=exon_coord[k,1]
        new_col3=exon_coord[k,2]
        col2=c(col2,new_col2)
        col3=c(col3,new_col3)
      }
      if(diff1==3)
      {
        if(k<3) temp=substring(amino_acid,(length(aa)+1),(length(aa)+length(new_col2)))
        else temp=substring(amino_acid,(length(aa)+3-k),(length(aa)+length(new_col2)+2-k))
        aa=c(aa,as.character(strsplit(temp,"")[[1]]))
      }
      else
      {
        if(k<3) temp=substring(amino_acid,(length(aa)),(length(aa)+length(new_col2)-1))
        else temp=substring(amino_acid,(length(aa)+2-k),(length(aa)+length(new_col2)+1-k))
        aa=c(aa,as.character(strsplit(temp,"")[[1]]))
      }
    }
    bedgraph=rbind(bedgraph,cbind(chromosome,col2,col3,"+",aa,geneID))
  }
  if(length(col2)!=length(col3) | length(col2)!=length(aa))
  {
    warning_geneID=c(warning_geneID,geneID)
  }
}
write.table(bedgraph,"cerevisiae_proteome.bedgraph",col.names=F, row.names=F, quote=F, sep="\t")
write.table(warning_geneID,"warning.txt",col.names=F, row.names=F, quote=F, sep="\t")
