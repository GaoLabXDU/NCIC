clear 
clc
cnv_data=importdata('data/KIRC_cnv.csv');
cnv_Matrix_data=importdata('data/KIRC_cnv_Matrix.csv');
cnv=cnv_data.data;
cnv_Matrix=cnv_Matrix_data.data;

m_data=importdata('data/KIRC_m.csv');
m_Matrix_data=importdata('data/KIRC_m_Matrix.csv');
m=m_data.data;
m_Matrix=m_Matrix_data.data;

me_data=importdata('data/KIRC_me.csv');
me_Matrix_data=importdata('data/KIRC_me_Matrix.csv');
me=me_data.data;
me_Matrix=me_Matrix_data.data;

d=0.25;
rank_cnv=RankGene(cnv,cnv_Matrix,d);
rank_m=RankGene(m,m_Matrix,d);
rank_me=RankGene(me,me_Matrix,d);

for i=2:8
	label=NCNMTF_3_type(cnv,m,me,rank_cnv,rank_m,rank_me,i);
	write_label=sprintf('%s%d%s','res/KIRC_cnv_m_me_',i,'.csv');
	dlmwrite(write_label,label);
end