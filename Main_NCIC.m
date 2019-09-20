clear 
clc

%input data contain expression data and Network data
%	****_cnv.csv is the expression data      	row is gene col is sample
%	****_cnv_Matrix.csv is the network data    	row and col are gene 

cnv_data=importdata('data/****_cnv.csv');
cnv_Matrix_data=importdata('data/****_cnv_Matrix.csv');
cnv=cnv_data.data;
cnv_Matrix=cnv_Matrix_data.data;

m_data=importdata('data/****_m.csv');
m_Matrix_data=importdata('data/****_m_Matrix.csv');
m=m_data.data;
m_Matrix=m_Matrix_data.data;

me_data=importdata('data/****_me.csv');
me_Matrix_data=importdata('data/****_me_Matrix.csv');
me=me_data.data;
me_Matrix=me_Matrix_data.data;

%input the parameter
d=0.25;

%each gene get different rank score 
rank_cnv=RankGene(cnv,cnv_Matrix,d);
rank_m=RankGene(m,m_Matrix,d);
rank_me=RankGene(me,me_Matrix,d);


%get the label for each K from 2~8
for i=2:8
	label=NCNMTF_3_type(cnv,m,me,rank_cnv,rank_m,rank_me,i);
	write_label=sprintf('%s%d%s','res/****_cnv_m_me_',i,'.csv');
	dlmwrite(write_label,label);
end
clear i label write_label

% if you have one type data or two type data,you can get the label like:
% for i=2:8
%	label=NCNMTF_1_type(cnv,rank_cnv,i);
%	write_label=sprintf('%s%d%s','res/****_cnv_',i,'.csv');
%	dlmwrite(write_label,label);
%end
%clear i label write_label
%
%
%for i=2:8
%	label=NCNMTF_1_type(cnv,m,rank_cnv,rank_m,i);
%	write_label=sprintf('%s%d%s','res/****_cnv_m_',i,'.csv');
%	dlmwrite(write_label,label);
%end
%clear i label write_label


