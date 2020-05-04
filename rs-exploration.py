'''
What will I be doing in this file?
- create new text file containinng fasta header and sequences separated by tabs
- entire sequence needs to be in a single line
- output needs to be gzipped
'''

'''
import sys

import gzip
from mimetypes import guess_type
from functools import partial

from Bio import SeqIO

input_file = sys.argv[1]
filename = str(input_file.split('.')[0])
print ("The file is " + filename) # prints file name (without extensions)
# this code is for opening file, came from: https://stackoverflow.com/questions/42757283/seqio-parse-on-a-fasta-gz
encoding = guess_type(input_file)[1]  # uses file extension
if encoding is None:
    _open = open
elif encoding == 'gzip':
    _open = partial(gzip.open, mode='rt')
else:
    raise ValueError('Unknown file encoding: "{}"'.format(encoding))


# this creates a sequence record
with _open(input_file) as f:
    
    for record in SeqIO.parse(f, 'fasta'):
        sys.stdout=open(filename + ".txt", "w")
        print(record)
        sys.stdout.close()
    
    for lines in f:
    	if '>' in lines:
    		header = lines
    	if '>' not in lines:
    		seq = lines #''.join(lines.replace('\n', ''))
    		print(header.replace('\n', '') + ',' + seq) 
    		# need to do something like this: https://www.geeksforgeeks.org/given-two-strings-find-first-string-subsequence-second/
    		# if string before the , matches, append sequences after ,
    		#for news in seq:
    			#sys.stdout=open(filename + ".csv", "w")
				#print(header.replace('\n', '') + ',' + seq)
        		#sys.stdout.close()


    	#seq = lines[::-1].replace("\n","",(lines.split('>')).count("\n")-1)[::-1]
    	#print(seq)

#close(input_file)
'''
