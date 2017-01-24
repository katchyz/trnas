### split tRNAs
### write to separate files for each anticodon
## K A Chyzynska, Jan 2017 ##

import string
import textwrap


f = open("/Users/kasia/Documents/PhD/scripts/tRNAs/medaka/all.fa", 'r')

## get headers
headers = []
header = ''
for line in f:
	line = line.strip()
	if line.startswith('>'):
		if header and seq:
			# save header and seq
			headers.append(header)
			print header
			print seq
		header = line.split()[4]+"_"+line.split()[5][1:4]
		seq = ''
	else:
		seq += line
	# for the last seq
	# save header and seq
	headers.append(header)


#######
## open files
headers = list(set(headers))

class A(object): pass

h = A()

for i in range(len(headers)):
	print headers[i]
	setattr(h, headers[i], open(headers[i]+'.txt', 'w'))


##########
## save fasta
f = open("/Users/kasia/Documents/PhD/scripts/tRNAs/medaka/all.fa", 'r')

header = ''
seq = ''
for line in f:
	line = line.strip()
	if line.startswith('>'):
		if header and seq:
			# remove introns (if any)
			seq = seq.translate(None,string.ascii_lowercase)
			# save header and seq
			getattr(h, header).writelines(['>', header, '_', chrom, '\n', textwrap.fill(seq, 60), '\n'])
		header = line.split()[4]+"_"+line.split()[5][1:4]
		print line.split()[2]
		chrom = line.split()[2]
		seq = ''
	else:
		seq += line


# for the last seq (skip on 1st line)
seq = seq.translate(None,string.ascii_lowercase)
getattr(h, header).writelines(['>', header, '_', chrom, '\n', textwrap.fill(seq, 60) , '\n'])



##########
## close handles
for a in dir(h):
	if not a.startswith('__'):
		print(a)
		getattr(h, a).close()






