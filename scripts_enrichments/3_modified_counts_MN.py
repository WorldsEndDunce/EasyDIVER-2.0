# Script to generate modified count files (one perselection)
# valid for Matilda Newton's project
# Requires bash file modified_counts_bash_MN.sh

# Use as: modified_counts_bash_MN.sh selection_in selection_neg selection_out output_file

# ATP:
# bash ../modified_counts_bash_MN.sh 3-5R5-Fus-in_S3_counts_aa_L.txt 4-5R5-Neg_S4_counts_aa_L.txt 5-5R5-ATP_S5_counts_aa_L.txt modified2/5R5-ATP_counts_aa_L_enr.txt
# bash ../modified_counts_bash_MN.sh 10-9R7-Fus-in_S10_counts_aa_L.txt 11-9R7-Neg_S11_counts_aa_L.txt 12-9R7-ATP_S12_counts_aa_L.txt modified2/9R7-ATP_counts_aa_L_enr.txt
# bash ../modified_counts_bash_MN.sh 17-16R5-Fus-in_S17_counts_aa_L.txt 18-16R5-Neg_S18_counts_aa_L.txt 19-16R5-ATP_S19_counts_aa_L.txt modified2/16R5-ATP_counts_aa_L_enr.txt 
# bash ../modified_counts_bash_MN.sh 24-20R5-Fus-in_S24_counts_aa_L.txt 25-20R5-Neg_S25_counts_aa_L.txt 26-20R5-ATP_S26_counts_aa_L.txt modified2/20R5-ATP_counts_aa_L.txt
# 
# GTP:
# bash ../modified_counts_bash_MN.sh 3-5R5-Fus-in_S3_counts_aa_L.txt 4-5R5-Neg_S4_counts_aa_L.txt 6-5R5-GTP_S6_counts_aa_L.txt modified2/5R5-GTP_counts_aa_L_enr.txt
# bash ../modified_counts_bash_MN.sh 10-9R7-Fus-in_S10_counts_aa_L.txt 11-9R7-Neg_S11_counts_aa_L.txt 13-9R7-GTP_S13_counts_aa_L.txt modified2/9R7-GTP_counts_aa_L_enr.txt
# bash ../modified_counts_bash_MN.sh 17-16R5-Fus-in_S17_counts_aa_L.txt 18-16R5-Neg_S18_counts_aa_L.txt 20-16R5-GTP_S20_counts_aa_L.txt modified2/16R5-GTP_counts_aa_L_enr.txt 
# bash ../modified_counts_bash_MN.sh 24-20R5-Fus-in_S24_counts_aa_L.txt 25-20R5-Neg_S25_counts_aa_L.txt 27-20R5-GTP_S27_counts_aa_L.txt modified2/20R5-GTP_counts_aa_L_enr.txt
# 
# ATPGTP:
# bash ../modified_counts_bash_MN.sh 3-5R5-Fus-in_S3_counts_aa_L.txt 4-5R5-Neg_S4_counts_aa_L.txt 7-5R5-ATPGTP_S7_counts_aa_L.txt modified2/5R5-ATPGTP_counts_aa_L_enr.txt
# bash ../modified_counts_bash_MN.sh 10-9R7-Fus-in_S10_counts_aa_L.txt 11-9R7-Neg_S11_counts_aa_L.txt 14-9R7-ATPGTP_S14_counts_aa_L.txt modified2/9R7-ATPGTP_counts_aa_L_enr.txt
# bash ../modified_counts_bash_MN.sh 17-16R5-Fus-in_S17_counts_aa_L.txt 18-16R5-Neg_S18_counts_aa_L.txt 21-16R5-ATPGTP_S21_counts_aa_L.txt modified2/16R5-ATPGTP_counts_aa_L_enr.txt 
# bash ../modified_counts_bash_MN.sh 24-20R5-Fus-in_S24_counts_aa_L.txt 25-20R5-Neg_S25_counts_aa_L.txt 28-20R5-ATPGTP_S28_counts_aa_L.txt modified2/20R5-ATPGTP_counts_aa_L_enr.txt

from time import time
import sys
from operator import itemgetter
import os


start = time()

if not os.path.exists("modified/"):
    os.makedirs("modified/")

print(sys.argv[1])
print(sys.argv[2])

fus_file  = sys.argv[1] # '3-5R5-Fus-in_S3_counts_aa_L.txt'  
neg_file = sys.argv[2] 		# '4-5R5-Neg_S4_counts_aa_L.txt' 
postsel_file = sys.argv[3] # '5-5R5-ATP_S5_counts_aa_L.txt' OR '6-5R5-GTP_S6_counts_aa_L.txt' OR '7-5R5-ATPGTP_S7_counts_aa_L.txt'

naming_ref   = sys.argv[5] # '7-5R5-ATPGTP_S7_counts_aa_L.txt'

out_file      = sys.argv[4] # '5R5-ATP_counts_aa_L_enr.txt' 

alph = sys.argv[6]

temp_file = 'modified/test.txt'

sel = postsel_file.split("-")[1]


files = [fus_file, neg_file, postsel_file, naming_ref]

illegal_5 = ['R', 'N', 'C', 'E', 'Q', 'H', 'I', 'L', 'K', 'F', 'S', 'T', 'W', 'Y']
illegal_9 = ['R', 'N', 'C', 'Q', 'H', 'I', 'K', 'F', 'W', 'Y']
illegal_16 = ['F', 'W', 'Y']

def n_illegal_5(string):
	count = 0 
	for i in illegal_5:
		count += string.count(i)
	return count

def n_illegal_9(string):
	count = 0 
	for i in illegal_9:
		count += string.count(i)
	return count

def n_illegal_16(string):
	count = 0 
	for i in illegal_16:
		count += string.count(i)
	return count
	
exc_9 = ['E', 'L', 'S', 'T']
exc_16 = ['R', 'N', 'C', 'Q', 'H', 'I', 'K']
exc_20 = ['F', 'W','Y']


all_dict = []
totals = []
uniques = []
max_len = 0
rnk = 0
for in_file in files:
	seqs   = []
	abunds = []
	degs = []
	fracs  = []
	ranks = []
	rnk = 0
	with open(in_file) as in_data:
		unique_line = next(in_data)
		total_line  = next(in_data)
		unique = int(unique_line.split('=')[1])
		total  = int(total_line.split('=')[1])
		totals.extend([total])
		uniques.extend([unique])
		next(in_data)
		for line in in_data:
			rnk += 1
			#print line
			seq = line.split()[0]
			abund = int(line.split()[1])
			deg = int(line.split()[3])
			frac  = abund / float(total)
			seqs.extend([seq])
			abunds.extend([abund])
			degs.extend([deg])
			fracs.extend([frac])
			ranks.extend([rnk])
			if len(seq) > max_len:
				max_len = len(seq)
		seqfit_list = [[seqs[i], ( abunds[i],fracs[i], ranks[i], degs[i])] for i in range(len(seqs))]
		seqfit_dict = dict(seqfit_list)
		all_dict.append(seqfit_dict)

print "max:", max_len

#print all_dict

out = open(temp_file,'w')
	
print >> out, str('number of unique sequences fus / neg / post = ') + str(uniques[0]) + str(" / ")+ str(uniques[1])  + str(" / ")+ str(uniques[2])
print >> out, str('total number of molecules fus / neg / post = ') + str(totals[0]) + str(" / ")+ str(totals[1])  + str(" / ")+ str(totals[2])
print >> out, str('rank').ljust(10),
print >> out, str('seq').ljust(max_len),
 
print >> out, str('a_cDNA').ljust(10), str('f_cDNA').ljust(20), 

print >> out, str('a_cDNA_l').ljust(10), str('f_cDNA_l').ljust(20), 
print >> out, str('a_cDNA_u').ljust(10), str('f_cDNA_u').ljust(20),  
 
print >> out, str(str('a_')+str(postsel_file.split("_")[0]).split("-")[-1]).ljust(10),
print >> out, str(str('f_')+str(postsel_file.split("_")[0]).split("-")[-1]).ljust(20),

print >> out, str(str('a_')+str(neg_file.split("_")[0]).split("-")[-1]).ljust(10),
print >> out, str(str('f_')+str(neg_file.split("_")[0]).split("-")[-1]).ljust(20),

print >> out, str(str('a_')+str(neg_file.split("_")[0]).split("-")[-1]+str("_l")).ljust(10),
print >> out, str(str('f_')+str(neg_file.split("_")[0]).split("-")[-1]+str("_l")).ljust(20),

print >> out, str(str('a_')+str(neg_file.split("_")[0]).split("-")[-1]+str("_u")).ljust(10),
print >> out, str(str('f_')+str(neg_file.split("_")[0]).split("-")[-1]+str("_u")).ljust(20),

print >> out, str(str('e_')+str(postsel_file.split("_")[0]).split("-")[-1]).ljust(20),
print >> out, str(str('e_')+str(postsel_file.split("_")[0]).split("-")[-1]+str("_l")).ljust(20),
print >> out, str(str('e_')+str(postsel_file.split("_")[0]).split("-")[-1]+str("_u")).ljust(20),


print >> out, str(str('e_neg')).ljust(20),
print >> out, str(str('e_neg')+str("_l")).ljust(20),
print >> out, str(str('e_neg')+str("_u")).ljust(20),

print >> out, str(str(str('e_')+str(postsel_file.split("_")[0]).split("-")[-1]) + ("/") + str(str('e_neg'))).ljust(20),

print >> out, str(str(str('e_')+str(postsel_file.split("_")[0]).split("-")[-1]) + ("/") + str(str('e_neg_u'))).ljust(20),
print >> out, str(str(str('e_')+str(postsel_file.split("_")[0]).split("-")[-1]) + ("/") + str(str('e_neg_l'))).ljust(20),


print >> out, str(str('ill')).ljust(10),
print >> out, str(str('n_ill')).ljust(10),
print >> out, str(str('from_red')).ljust(20),

print >> out, str(str('deg')).ljust(10)

new_total = 0
new_unique = 0

for seq in all_dict[2]:
	c_post = all_dict[2][seq][0]
	f_post = all_dict[2][seq][1]
	deg_post = all_dict[2][seq][3]
	
	try:
		rank_post = all_dict[3][seq][2]
	
	except KeyError:
		rank_post = str("-")
	
	try:
		c_fus  = all_dict[0][seq][0]
		f_fus  = all_dict[0][seq][1]
		
		# numbers from low counts excel file
		if c_fus == 1:
			c_fus_l = 1
			c_fus_u = 3
			f_fus_l = float(1)/totals[0]
			f_fus_u = float(3)/totals[0]
		if c_fus == 2:
			c_fus_l = 1
			c_fus_u = 5
			f_fus_l = float(1)/totals[0]
			f_fus_u = float(5)/totals[0]
		if c_fus == 3:
			c_fus_l = 1
			c_fus_u = 7
			f_fus_l = float(1)/totals[0]
			f_fus_u = float(7)/totals[0]
		if c_fus == 4:
			c_fus_l = 1
			c_fus_u = 8
			f_fus_l = float(1)/totals[0]
			f_fus_u = float(8)/totals[0]
		if c_fus == 5:
			c_fus_l = 1
			c_fus_u = 10
			f_fus_l = float(1)/totals[0]
			f_fus_u = float(10)/totals[0]
		if c_fus == 6:
			c_fus_l = 2
			c_fus_u = 11
			f_fus_l = float(2)/totals[0]
			f_fus_u = float(11)/totals[0]
		if c_fus == 7:
			c_fus_l = 2
			c_fus_u = 13
			f_fus_l = float(2)/totals[0]
			f_fus_u = float(13)/totals[0]
		if c_fus == 8:
			c_fus_l = 3
			c_fus_u = 14
			f_fus_l = float(3)/totals[0]
			f_fus_u = float(14)/totals[0]
		if c_fus == 9:
			c_fus_l = 4
			c_fus_u = 15
			f_fus_l = float(4)/totals[0]
			f_fus_u = float(15)/totals[0]
		if c_fus == 10:
			c_fus_l = 4
			c_fus_u = 17
			f_fus_l = float(4)/totals[0]
			f_fus_u = float(17)/totals[0]
		if c_fus > 10:
			c_fus_l = str("-")
			c_fus_u = str("-")
			f_fus_l = str("-")
			f_fus_u = str("-")

	except KeyError:
	
		c_fus = 1
		f_fus = float(1)/totals[0]  # Change this line to modify fraction used if not found
		
		c_fus_l = str("-")
		c_fus_u = str("-")
		f_fus_l = str("-")
		f_fus_u = str("-")



	try:
		c_neg  = all_dict[1][seq][0]
		f_neg  = all_dict[1][seq][1]
		
		if c_neg == 1:
			c_neg_l = 1
			c_neg_u = 3
			f_neg_l = float(1)/totals[0]
			f_neg_u = float(3)/totals[0]
		if c_neg == 2:
			c_neg_l = 1
			c_neg_u = 5
			f_neg_l = float(1)/totals[0]
			f_neg_u = float(5)/totals[0]
		if c_neg == 3:
			c_neg_l = 1
			c_neg_u = 7
			f_neg_l = float(1)/totals[0]
			f_neg_u = float(7)/totals[0]
		if c_neg == 4:
			c_neg_l = 1
			c_neg_u = 8
			f_neg_l = float(1)/totals[0]
			f_neg_u = float(8)/totals[0]
		if c_neg == 5:
			c_neg_l = 1
			c_neg_u = 10
			f_neg_l = float(1)/totals[0]
			f_neg_u = float(10)/totals[0]
		if c_neg == 6:
			c_neg_l = 2
			c_neg_u = 11
			f_neg_l = float(2)/totals[0]
			f_neg_u = float(11)/totals[0]
		if c_neg == 7:
			c_neg_l = 2
			c_neg_u = 13
			f_neg_l = float(2)/totals[0]
			f_neg_u = float(13)/totals[0]
		if c_neg == 8:
			c_neg_l = 3
			c_neg_u = 14
			f_neg_l = float(3)/totals[0]
			f_neg_u = float(14)/totals[0]
		if c_neg == 9:
			c_neg_l = 4
			c_neg_u = 15
			f_neg_l = float(4)/totals[0]
			f_neg_u = float(15)/totals[0]
		if c_neg == 10:
			if alph == '20aa':
				c_neg_l = 4
				c_neg_u = 16
				f_neg_l = float(4)/totals[0]
				f_neg_u = float(16)/totals[0]
			else:
				c_neg_l = 4
				c_neg_u = 17
				f_neg_l = float(4)/totals[0]
				f_neg_u = float(17)/totals[0]
		if c_neg > 10:
			c_neg_l = str("-")
			c_neg_u = str("-")
			f_neg_l = str("-")
			f_neg_u = str("-")
		
	except KeyError:
		c_neg = 1
		f_neg = float(1)/totals[1]  # Change this line to modify fraction used if not found
	
		c_neg_l = str("-")
		c_neg_u = str("-")
		f_neg_l = str("-")
		f_neg_u = str("-")
		
	# Calculate and adjust enrichment in positive and negative pools
	enr_post   = f_post / f_fus     # Enrichment due to selection
	enr_neg  = f_neg / f_fus     # Enrichment due to selection
	
	if type(f_fus_l) == float:
	
		enr_post_l   = f_post / f_fus_u     # Enrichment due to selection - upper
		enr_post_u   = f_post / f_fus_l     # Enrichment due to selection - lower

		if type(f_neg_l) == float:
		
			enr_neg_l  = f_neg_l / f_fus_u     # Enrichment due to selection - upper
			enr_neg_u  = f_neg_u / f_fus_l     # Enrichment due to selection - lower
		
		else:
		
			enr_neg_l  = f_neg / f_fus_u     # Enrichment due to selection - upper
			enr_neg_u  = f_neg / f_fus_l     # Enrichment due to selection - lower
	else:
		
		if type(f_neg_l) == float:
		
			enr_neg_l  = f_neg_l / f_fus     # Enrichment due to selection - upper
			enr_neg_u  = f_neg_u / f_fus     # Enrichment due to selection - lower
		
		else:
		
			enr_neg_l  = str("-")     # Enrichment due to selection - upper
			enr_neg_u  = str("-")     # Enrichment due to selection - lower

	
		enr_post_l   = str("-")     # Enrichment due to selection - upper
		enr_post_u   = str("-")     # Enrichment due to selection - lower
		
		#enr_neg_l  = str("-")     # Enrichment due to selection - upper
		#enr_neg_u  = str("-")     # Enrichment due to selection - lower
	
	if alph == '5aa':
	
		#illegal:
		if any((res in illegal_5) for res in seq[0:-3]):
			seq_ill = 'yes'
		else:
			seq_ill = 'no'
		if n_illegal_5(seq[0:-3]) > 0:
			seq_n_ill = n_illegal_5(seq[0:-3])
		else:
			seq_n_ill = str('-')
		# from reduced alphabet:
		seq_from_red = str('-')
	
	
	if alph == '9aa':
		
		#illegal:
		if any((res in illegal_9) for res in seq[0:-3]):
			seq_ill = 'yes'
			seq_from_red = str('no')
		else:
			seq_ill = 'no'
				
			# from reduced alphabet:
			if any((res in exc_9) for res in seq[0:-3]):
				seq_from_red = str('no')
			else:
				seq_from_red = str('from_5AA')
				
		if n_illegal_9(seq[0:-3]) >0:
			seq_n_ill = n_illegal_9(seq[0:-3])
		else:
			seq_n_ill = str('-')
	
	if alph == '16aa':
		
		#illegal:
		if any((res in illegal_16) for res in seq[0:-3]):
			seq_ill = 'yes'
			seq_from_red = str('no')
		else:
			seq_ill = 'no'
		
			# from reduced alphabet:
			if any((res in exc_16) for res in seq[0:-3]):
				seq_from_red = str('no')
			else:
				if any((res in exc_9) for res in seq[0:-3]):
					seq_from_red = str('from_9AA')
				else:
					seq_from_red = str('from_5AA')
					
		if n_illegal_16(seq[0:-3]) >0:
			seq_n_ill = n_illegal_16(seq[0:-3])
		else:
			seq_n_ill = str('-')
			
	if alph == '20aa':
		#illegal:
		seq_ill = '-'
		seq_n_ill = str('-')
		# from reduced alphabet:
		if any((res in exc_20) for res in seq[0:-3]):
			seq_from_red = str('no')
		else:
			if any((res in exc_16) for res in seq[0:-3]):
				seq_from_red = str('from_16AA')
			else:
				if any((res in exc_9) for res in seq[0:-3]):
					seq_from_red = str('from_9AA')
				else:
					seq_from_red = str('from_5AA')


	# Write data to file
	#out.write('%s\t%.3e\t%d\t%.3e\t%d\t%.3e\t%d\t%.3f\n' % (seq,f_fus,c_fus,f_post,c_post,enr))
	
	#out.write('%s\t%.3e\t%d\t%.3e\t%d\t%.3f\n' % (seq,f_fus,c_fus,f_post,c_post,enr))
	#if float(f_post) > 0.0000711526: 				# corresponding f_L for 10 count reads for 20AA-ATP (which is f_L_cutoff = 10 / 140543 = 0.0000711526)
	if float(f_post) > 0.0000711526:
		#print "f_post", f_post
		new_total += int(c_post)
		new_unique += 1	
		
		##### name 
		 
		print >> out, str(str(sel)+str("-")+str("V")+str(rank_post)).ljust(10),
		
		##### sequence
		
		print >> out, str(seq).ljust(max_len), 
		
		##### abd anf freq fus
		
		print >> out, str(c_fus).ljust(10), 
		print >> out, str('%.10f' % float(f_fus)).ljust(20),
		
		##### abd and freq fusion LOWER and UPPER
		
		print >> out, str(c_fus_l).ljust(10), 
		if type(f_fus_l) == float:
			print >> out, str('%.10f' % float(f_fus_l)).ljust(20), 
		else:
			print >> out, str(f_fus_l).ljust(20), 
		
		print >> out, str(c_fus_u).ljust(10), 
		if type(f_fus_u) == float:
			print >> out, str('%.10f' % float(f_fus_u)).ljust(20),
		else:
			print >> out, str(f_fus_u).ljust(20),
		
		##### abd and freq post
		
		print >> out, str(c_post).ljust(10), 
		print >> out, str('%.10f' % float(f_post)).ljust(20), 
		
		##### abd and freq neg
		
		print >> out, str(c_neg).ljust(10), 
		print >> out, str('%.10f' % float(f_neg)).ljust(20), 
		
########################################################################### new
		##### abd and freq neg LOWER and UPPER
		
		print >> out, str(c_neg_l).ljust(10), 
		if type(f_neg_l) == float:
			print >> out, str('%.10f' % float(f_neg_l)).ljust(20), 
		else:
			print >> out, str(f_neg_l).ljust(20), 
		
		print >> out, str(c_neg_u).ljust(10), 
		if type(f_neg_u) == float:
			print >> out, str('%.10f' % float(f_neg_u)).ljust(20),
		else:
			print >> out, str(f_neg_u).ljust(20),
###########################################################################		
		
		##### enr post
		
		print >> out, str('%.10f' % float(enr_post)).ljust(20),
		
		##### enr post LOWER and UPPER
		
		if type(enr_post_u) == float:
			print >> out, str('%.10f' % float(enr_post_l)).ljust(20),
			print >> out, str('%.10f' % float(enr_post_u)).ljust(20),
		else:
			print >> out, str(enr_post_l).ljust(20),
			print >> out, str(enr_post_u).ljust(20),
		
		##### enr neg
		
		print >> out, str('%.10f' % float(enr_neg)).ljust(20),


		##### enr neg LOWER and UPPER ####################################### just modified	

		
		if type(enr_neg_u) == float:
			print >> out, str('%.10f' % float(enr_neg_l)).ljust(20),
			print >> out, str('%.10f' % float(enr_neg_u)).ljust(20),
		else:
			print >> out, str(enr_neg_l).ljust(20),
			print >> out, str(enr_neg_u).ljust(20),
		###########################################################################		

		##### rel enr
		
		if float(enr_neg) > 0:
			print >> out, str('%.10f' % float(enr_post/enr_neg)).ljust(20),
		else:
			print >> out, str("-").ljust(20),

###########################################################################		
			
		# rel enr LOWER and UPPER
		
		if type(enr_neg_u) == float:
			print >> out, str('%.10f' % float(enr_post/enr_neg_u)).ljust(20),
		else:
			print >> out, str("-").ljust(20),
		if type(enr_neg_u) == float:
			print >> out, str('%.10f' % float(enr_post/enr_neg_l)).ljust(20),
		else:
			print >> out, str("-").ljust(20),

###########################################################################		
		
		# illegals and red alphabet
		
		print >> out, str(seq_ill).ljust(10),
		print >> out, str(seq_n_ill).ljust(10),
		print >> out, str(seq_from_red).ljust(10),
		
		print >> out, str(deg_post).ljust(10)
		
print (time()-start)

out.close()


out = open(temp_file,'r')
out2 = open(out_file,'w')
count_line = 0
for line in out:
	count_line += 1
	if count_line == 1:
		print >> out2, str(line.split("\n")[0] + str( " | number of unique sequences post = ") + str(new_unique))
	if count_line == 2:
		print >> out2, str(line.split("\n")[0] + str( " | total number of molecules = ") + str(new_total))
	if count_line != 1 and count_line != 2:
		print >> out2, line,
out.close()
out2.close()

os.remove(temp_file)
