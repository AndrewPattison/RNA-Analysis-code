

###This functon asks for a number input###

def getinput():
	num = raw_input("Please enter the number of tubes you will be using in your assay:")
	return num

#Set a variable to the result of our input
num = getinput() 

### Checks that the input can be made an int.###
def RepresentsInt(s):
    try: 
        int(s)
        return True
    except ValueError:
        return False
#If input is not an int, runs getinput again.
while RepresentsInt(num) ==False:
	print '\nThat was not a positive intiger!\n'
	num = getinput() 
else: 
	num = num

#Warning if intiger entered is not positive
if int(num) <= 0:
	print 'Warning: non positive int used, returned values will be negatives'





###Asks for an assay selection###	
def get_assay():
	assay = input('Which assay are you making the mix for? \n 1. EM-Pat \n 2. T12VN-PAT size control cDNA \n 3. EM-PAT PCR \n ')	
	return assay

#Sets selction to a variable
selection = get_assay()

#List of acceptable selections
answers = [1,2,3]

#Loop to make sure selection is one of approved selections
while selection not in answers:
	print '\nPlease type either 1, 2 or 3\n'
	selection = get_assay()


list_for_EM_PAT = [['RNA/dH20', 7], ['PAT assay primer', 1], ['dH2O', 4], ['5X SSIII buffer', 4], ['0.1M DTT', 1], ['dNTPs (10mM each)', 1], ['RnaseOut', 1], ['Klenow polymerase', 1]]
list_for_T12VN_PAT = [['RNA/dH20', 7], ['T12VN PAT size control oligo (100 uM)', 1], ['dH2O', 4], ['5 x superscript reaction buffer; (Invitrogen)', 4], ['0.1M DTT', 1], ['dNTPs (10mM each)', 1], ['RnaseOut', 1]]
list_for_EM_PAT_PCR = [['Amplitaq Gold 360 Master Mix', 10], ['PAT-Assay oligos', 0.2], ['Gene of interest oligos', 0.2], ['diluted cDNA (EM-PAT & T12VN-PAT)', 5], ['dH2O', 4.6]]

	
if selection == 1:
	print 'Volumes of solutions to mixed for EM PAT: \n \n '
	total = 0
	for pair in list_for_EM_PAT:
		total += (pair[1]*0.5) + (pair[1]* float(num))
		print  '%s:' %(pair[0]), '%s' % (float (num)*float (pair[1]) + 0.5*float (pair[1]) ), 'ul\n'
	print 'Total volume of EM_PAT master mix is:', total, 'ul\n'


if selection == 2:
	print 'Volumes of solutions to mixed for T12VN PAT: \n \n '
	total = 0
	for pair in list_for_T12VN_PAT:
		total += (pair[1]*0.5) + (pair[1]* float(num))
		print  '%s:' %(pair[0]), '%s' % (float (num)*float (pair[1]) + 0.5*float (pair[1]) ), 'ul\n'
		
	print 'Total volume of T12VN PAT master mix is:', total, 'ul\n (for 19 ul of master mix per tube) 1ul of Superscript III\n added (per tube) later\n'


if selection == 3:
	print 'Volumes of solutions to mixed for EM PAT PCR: \n \n '
	total = 0
	for pair in list_for_EM_PAT_PCR:
		total += (pair[1]*0.5) + (pair[1]* float(num))
		print  '%s:' %(pair[0]), '%s' % (float (num)*float (pair[1]) + 0.5*float (pair[1]) ), 'ul\n'
	print 'Total volume of EM PAT PCR master mix is:', total, 'ul\n'


#go_again = input('Would you like to do another assay?' '(y or n)')

#if go_again == 'y':
	 
	
#else:
	print 'goodbye' 