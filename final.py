import numpy as np
import time 
import random
import matplotlib.pyplot as plt 
from collections import Counter
import mimetypes


## string preparation functions: seq_generator, back to file .

# the sequence generator converts files to binary or generate a random sequence of a given length
def seq_generator(input,input_mode):

     # input mode 2 is a file import

    if int(input_mode) == 2:
        path = input
        try:
            with open(path, "rb") as file:
                w_in_bytes = file.read()
                w_in = []
                for byte in w_in_bytes:
                    bits = format(byte, '08b')
                    w_in.extend(int(bit) for bit in bits)
                print("File converted to binary successfully.")
        except FileNotFoundError:
            print("File not found. Please enter a valid file path.")
        except Exception as e:
            print(f"Error while processing the file: {e}")
        n = len(w_in)+2
        k = int(np.ceil(2*(1+np.log2(n))))
        if k % 2 != 0:    # to avoid length related issues
            k=k+1   
 
    # input mode 1 is a random sequence of a chosen length 
   
    elif int(input_mode) == 1:
        n = int(input)+2
        k = int(np.ceil(2*(1+np.log2(n))))
        if k % 2 != 0:    # to avoid length related issues
            k=k+1   
        w_in = []
        

       # seed_value = 42
       # random.seed(seed_value)
 #creating a random sequence. note by definiotion random would note yield many repetitions, especialy in short strings. 
        w_in = [random.randint(0,1) for _ in range(n-2)]


        # to create a higher rate of repetitions the following alteration can be made :
        for j in range(len(w_in)):
            if j % 4 == 0:
                w_in[j]=0

    return w_in,k


# the back to file function is optional for use in the encoding and decoding modes.
#in encoding it will autamaticlly be saved as txt file containing the encoded string, which can later be sent for synthesis .
#in decoding it will be saved in a f0rmat according to build in indicaters in the binary data.

def back_to_file(w,path):

    bits = w
    bytes_data = bytes([int(''.join(map(str, bits[i:i+8])), 2) for i in range(0, len(bits), 8)])
    with open(path, "wb") as file:
        file.write(bytes_data)
 


### encoding functions:
## 1)encoding support functions: pad marker, eliminate, check duplicates, generate binary words
## 2)enc1 - alg 0 :(dynamic Rabin Karp hash + hash)/ alg 1: (original + info)/ alg 2: (Boyer Moore + preprocessing 1 + preprocessing 2 + compute array )
## 3)enc2 , enc 3



## 1)

#pad marker function is the initiall stage of encoding, it adds the sufix and prefix markers.
def pad_marker(w,k):
    w_in=[]
    suf_marker = [1]
    suf_marker.extend([0]*int(np.ceil(k/2)))
    w_in.append(0)
    w_in.extend(w)
    w_in.extend(suf_marker)
    return w_in


# eliminate function is used to eliminate the repeted window as well as retaining the "adress" of both windsows 
# later in decoding the second adress will be used to obtain the eliminated sequence and the first as the location in which it should be reinserted 

def eliminate(w,idx1,idx2,k):
    # converting position number (decimal) to binary representation with length k/2 -1 (= log2(n))
    idx_1b=(bin(idx1)[2:]).zfill(int(np.ceil(k/2))-1)
    idx_1b =[int(x) for x in idx_1b]  
    idx_2b =(bin(idx2)[2:]).zfill(int(np.ceil(k/2))-1)
    idx_2b =[int(x) for x in idx_2b] 
    #making a new string with a marker,the indices and the original string with the first window removed
    w_temp = w[:idx1]+w[idx1+k:]    
    w=[0] 
    w.extend(idx_1b+idx_2b+w_temp)
    return w 


# check duplicate function is used to verify the operation of the encoder.
def check_duplicate(w, k):
    seen_lists = set()
    for i in range(len(w) - k + 1):
        window = w[i:i + k]
        window_tup = tuple(window)
        if window_tup in seen_lists:
            return True  # Found a duplicate window
        else:
            seen_lists.add(window_tup)
    return False  # No duplicate windows were found

# the generate_binary_word function is used during expansion (enc3), it calculates all possible adittions to the string 
#It employs a recursive construction of a binary tree, starting from a word of 0 bits and branching out by adding both '0' and '1' in each iteration for a total of k iterations. This process generates all possible permutations and combinations of k bits
def generate_binary_words(k):
   #
    def generate_words(current_word, remaining_length):
        if remaining_length == 0:   #
            binary_words.append(current_word.copy())  
            return
        # branch out '0'- branch :
        current_word.append(0)                                #current word --> current word+'0'
        generate_words(current_word, remaining_length - 1)    #k --> k-1
        current_word.pop()           
        #branch out '1' - branch: 
        current_word.append(1)                                 #current word --> current word+'1'
        generate_words(current_word, remaining_length - 1)     #k -->k-1
        current_word.pop()  
    
    binary_words = []
    generate_words([], k) # initial activation
    return binary_words

## 2)

# Rabin Karp hash is the decimal value of the binary stirng. it can be easily calculated given the hash of the previous window whos k-1 suffix is identical to the current windows k-1 prefix.
def R_K_hash (msb,lsb,last_hash,k):
    new_hash = 2*(last_hash-2**(k-1)*msb)+lsb
    return new_hash 

#alg_0: dynamic search using Rabin Karp hash. - default.
# dynamically creating a dictionary of new windows, the hash beeing the key of each window and the value is the starting position(index) of each window 
# when a key is identical to a key that already appears in the dictionary a repetitioin is identified
def r_k_enc1_dynamic(w,k):
    reps = 0
    flag=1 # lest run yielded a repeat , if flag is set to 0 and not changed (no repetitions found) then the string has not been changed and the proccess is over.
    while flag==1:
        flag = 0
        seen_windows = {}
        first_idx = -1  
        second_idx = -1

        for j in range(len(w) - k+ 1):
            if j==0:
                # first entery's hash must be calculated independetly 
                hash_j= ''.join(str(e) for e in w[0:k])
                hash_j = int(hash_j,2)
            else:  
                # following hash can be calculated in a less computationally complex way using R_K_hash function
                hash_j = R_K_hash(w[j-1],w[j+k-1],hash_j,k) 
            if hash_j in seen_windows.keys():
                # repetition found - updating the indices 
                if first_idx == -1 or (seen_windows[hash_j] < first_idx):  # the condition makes sure we only identify the first repetition i.e., the first window from the left that has another accurance later in the string.
                    first_idx = seen_windows[hash_j]  # first index takes the value of the index of the existing hash (privious appearance)  
                    second_idx = j                    # second index takes the value of the current index      
            else:
                seen_windows[hash_j] = j # new hash - saving the index as value with thenew hash as it's key.
       
        if first_idx >= 0: # two identical windows were found, as the index has been updated
            flag =1 # the string is beeing changed by the eliminate function therefor another run is needed
            reps = reps+1 #counts the number of iteration for later analysis 
            w = eliminate(w,first_idx,second_idx,k)
        if flag==0:
            break
    return w,reps
   


#alg_1 : the original search algorithm - preprocessing   
# The following search algorithm is my original creation
# for more information regarding the method of operation and the heuristics I used please see explanation in READ.ME.

# obtains information for original search algorithm's heurictics
def info(w,k):
    suffix_array =[]
    sum_indicator = []
    for j in range (len(w)-k+1):
        suffix_array.append(w[int(j+k/2):int(j+k)])
        sum_indicator.append(sum(w[int(j+k/2):int(j+k)]))
    return suffix_array,sum_indicator

 
def original_enc1(w,k):
    reps = 0
    lw=len(w)
    (suffix_array,sum_indicator)= info(w,k) 
    flag = 0
    i=0
    while i in range (len(w)-k):
        if flag==1:    # a repetition was found 
            w = eliminate(w,idx1,idx2,k)
            reps += 1
            i=0
            (suffix_array,sum_indicator)= info(w,k)
            flag=0
        else:
            # searching for windows identical to the pattern,
            # for each pattern w[i:i+k], the search starts at i+1 and scour a substring with all the windows that appears after the pattern in the string.
            suffix_match =[s for s in range(i+1,len(w)-k) if suffix_array[i]==suffix_array[s]] # checks for for heuristic (bad suffix)  
            if len(suffix_match)==0:
                i = i+int(k/2)+1 # skip k/2 + 1 patterns 
                continue
            sum_match = [j for j in suffix_match if sum_indicator[j] == sum_indicator[i]] # checks for second heuristic (sum indicator)
            if len (sum_match)==0:
                i=i+1 # skip to next window
                continue
            else:
                first_bit = [j for j in sum_match if w[i]==w[j]]
                # bitwise comparison of prefix only, as the suffix of the remaining windows has allready been compared and found to match the suffix of the pattern.
                prefix_len = int(np.ceil(k/2))
                repeats =[]
                for j in first_bit:
                    for p in range(1,prefix_len):
                        if w[i+p] != w[j+p]:
                            break  # Stop the inner loop if a mismatch is found
                    else:
                        repeats.append(j)

                if len(repeats)>0: 
                    flag =1
                    [idx1,idx2]= [i,repeats[0]]   #several repetitions of the pattern may exist, the first element in the repeats list is the first (from the left) appearance of the pattern  
                else:
                    i=i+1
    return w,reps


# The following search algorithm is a one-legged Boyer Moore.
#  The first of the two heuristics traditionally used in BM is the bad character rule which relies on having lots of different types of characters
#  as binary only has two it is not very useful in binary strings, therefore i only implemented the second heuristic, which is the good suffix rule.
# The good suffix rule shifts the pattern in comparison to the text such that a suffix t of the section in w being compared to the pattern is matched again or a suffix of it is matched with a prefix of the patterns 
# This is done using arrays for shifts and border positions (both suffix and prefix) which are calculated in the preprocessing functions.
# The preprocessing functions are taken from geeks for geeks at https://www.geeksforgeeks.org/boyer-moore-algorithm-good-suffix-heuristic/

# the good suffix heuristic is divided to two cases the first is a match of t to a substring in the pattern and the second is a match of a suffix of t to a prefix of p 



def case_2(shift, border_pos, pat, m):
    j = border_pos[0]   # the starting position of a substring that is both a suffix and a prefix of p (if there isnt one it is set to m+1)
    for i in range(m + 1):
        if shift[i] == 0:
            shift[i] = j
        if i == j:
            j = border_pos[j]

def case_1(shift, border_pos, pat, m):
    i = m
    j = m + 1
    border_pos[i] = j
    while i > 0:
        while j <= m and pat[i - 1] != pat[j - 1]:  # the charachter in i-1 is different than the one the mismeched in j-1 so when in j shift to i
            if shift[j] == 0:   # updating the shift
                shift[j] = j - i
            j = border_pos[j]   # the position of the first appearance of a substring that is both a prefix and a suffix of the suffix of p that starts in positon j
        i -= 1    # going backwards in the pattern 
        j -= 1
        border_pos[i] = j



def b_m_enc1(w,k):
    i=0 
    reps = 0 
    while i in range (len(w)-k):
        flag = 0
        pat =w[i:i+k]
        text = w[i+1:]
        s = 0   #intializing the shift
        m = len(pat)
        n = len(text)
        # initializing the shift and the border position arrays
        border_pos = [0] * (m + 1)
        shift = [0] * (m + 1)
        # calculating the shift and the border position arrays for the current string
        case_1(shift, border_pos, pat, m)
        case_2(shift, border_pos, pat, m)
        while s <= n - m:
            j = m - 1  # boyer moore compares from the last position of the pattern backwards
            # finding a suffix of p that matches a substring in the text.
            while j >= 0 and pat[j] == text[s + j]:
                j -= 1   #going backwards 
            if j < 0:   # the entire pattern has been matched - repeat found 
                idx2 = s+i+1
                w = eliminate(w,i,idx2,k)
                reps += 1
                i=0
                flag=1
                break
            else:
                s += shift[j + 1] # shifting the pattern to next appearance of t 
        if flag == 0:
            i=i+1
    return w,reps



## 3) 

# the second part of the encoder finds occurances of sufiix marker prior to suffix itself and eliminates them.
def enc_2(w,k):
    change_flag = False
    stop = 0
    while stop == 0:
        c=0
        for i in range(0,len(w)-int(k/2)): 
            if w[i:i+int(k/2)] == [0]*int(k/2):     # a suffix marker sequance found in the (i-1)'th position
                change_flag = True
                idx_suf = (bin(i)[2:]).zfill(int(k/2)-1)
                idx_suf =[int(x) for x in idx_suf] 
                w_temp = w[:i]+w[i+int(k/2):]
                w =[1]
                w.extend(idx_suf+w_temp)
                c=1
                continue  
        if c == 0 :
            stop = 1
    return w,change_flag

# third and final stage of encoding is the expansoin. this function adds bits to an encoded sequence so it will of fixed length(n). 
def enc_3(w,n,k):
    n_tag = len(w)
    diff = n-n_tag 
    addition_options =  generate_binary_words(diff) # calculates all the binary strings of the desired length
    num_opt = len(addition_options) 
    for i in range(num_opt):
        # each option for the final padded string is verified by the check duplicate function to have no repetitions, the first one that meets this demand will be the final output. 
        w_opt = w + addition_options[i]   
        if check_duplicate(w_opt,k)==False:  
            return w_opt
    raise Exception("w already contains all the binary words of length k.- expansion could not be preformed")

### decoding functions:
## remove add-ons, dec 1, dec 2

## remove add-ons funcion removes irelevent bits only retaining the sunstring of bits appearing prior to suffix marker. 
# it does so by differentiating between two states of the encoded string (prior to expansion in enc 3):
# The first state: encoded string was long enough(>n), enc3 was not activated. in this case suffix marker will be cut short and some of the '0''s at the end will be discarded
# The second state : encoded string was too short(<n), enc 3 was activated. in this case the suffix marker will appear as a whole but not at the end of the string  
def remove_addons(w_in_dec,n,k):
    marker_len = int(np.ceil(k/2))
    flag = 0
    position_zeros = n+1
    for i in range (1,n-marker_len+1):
        c=0
        if flag == 1:
            break
        for tau in range (0,marker_len):    #cheacking for the first appearance of suffix marker (case 2)
            if w_in_dec[i+tau] == 0:
                c+=1
            else:
                break
        if c == marker_len:            #  first suffix marker found . this is the point where the string that holds relevant data ends
            position_zeros = i 
            flag = 1   # whole suffix marker found we can exit the search 
    if flag == 0:     # whole suffix marker not found- encoded string did not go throgh expantion we must look for the begining of suf marker by looking for last'1' (case 1)
        j=1
        while j >=0:
            if w_in_dec[n-j] == 1:
                position_last1 = n-j    
                break
            else:
                j=j+1
        w = w_in_dec[0:position_last1-0]     # case 1 
    else:
        w = w_in_dec[0:position_zeros-1]      # case 2
    return w

## dec 1 reverses the operations made by enc 1  
def dec_1(w,k):
        #identifieng the adresses of the tuple
        idx_len = int(np.ceil(k/2))-1
        idx_1 = w[1:idx_len+1]
        idx_1 = ''.join(str(e) for e in idx_1)
        idx_1 = int(idx_1,2)
        idx_2 = w[idx_len+1:2*idx_len+1]
        idx_2 = ''.join(str(e) for e in idx_2)
        idx_2 = int(idx_2,2)

        # im case there was an overlap between the windows a part of the second window was deleted while deleting the first window
        # in this case the non ovelaping part is the period since we know the end of window1 is identical to the end of window2 and also overlaps with the begining of window 2       
        period = idx_2 -idx_1
        if period < k:
            u = w[idx_1+k-1:idx_1+k+period-1]
            # in case the end of the sequence was the suffix marker, which has already been removed
            if idx_2+k>len(w):
                u.extend([1]+[0]*(idx_2+k-len(w)-2))
            insertion_window = u*(int(k/period)+1)
            insertion_window =insertion_window[-k:] 
        else:
            insertion_window = w[idx_2-1:idx_2+k-1]
            if idx_2+k>len(w):
                insertion_window.extend([1]+[0]*(idx_2+k-len(w)-2))
      
        w_temp = w[k-1:] 
        w = w_temp[:idx_1]+ insertion_window +w_temp[idx_1:]
        return w

## dec 2 reverses the operations made by enc 2 
def dec_2(w,k):
        marker_len = int(np.ceil(k/2)+1)
        # identifieng the address of sufix marker sequence
        idx_suf = w[1:marker_len-1]
        idx = ''.join(str(e) for e in idx_suf)
        idx = int(idx,2)
        #reinstating the marker sequence
        w = w[marker_len-1:]
        insertion_window = [0]*(marker_len-1)
        w_out = w[:idx]+insertion_window+w[idx:] 
        return w_out



### the encoder decoder control functions:

## the encoder control function is in charge of preformaing the different operations of encoding (enc1- repetition elimination, enc2- suffix marker deduplication, enc3 - expansion ) in the following oreder:
# enc1 and enc2 must be preformed alternately until there are no apperances of both repetitions and suffix marker
# only then enc3 is employed if needed
def encoder(w_in,alg_num,k):
 
    time_line = []
    tot_reps = 0
    n = len(w_in)+2
  
    if alg_num ==0:
            func = r_k_enc1_dynamic
    else: 
        # in the test operating mode the user can choose a search algorithm. the output of all of them will be the same, however the run-time differs.
        func_dic=[r_k_enc1_dynamic,original_enc1,b_m_enc1]
        func = func_dic[int(alg_num)-1]
   
    w_padded = pad_marker(w_in,k)
   
    start_1 = time.time()
    w_1,reps= func(w_padded,k)
    time_line.append(['Enc 1',time.time()-start_1])
    tot_reps+=reps
   
    if reps == 0:  #input string has no k-repeats
       
        start_2 = time.time()
        (w_2,change_flag) = enc_2(w_1,k)   # even if the input string has no k-repeats the suffix marker must still be unique 
        time_line.append(['Enc 2',time.time()-start_2])
     
        if change_flag :
          
            start_1 = time.time()                   # changes to the string can result in a new k-repeat
            w_1,reps= func(w_2,k)
            time_line.append(['Enc 1',time.time()-start_1])
            tot_reps+=reps
    
    while reps > 0:                          # the alternating process decribed above
        
        start_2 = time.time()
        (w_2,change_flag) = enc_2(w_1,k)   
        time_line.append(['Enc 2',time.time()-start_2])
       
        if change_flag :
          
            start_1 = time.time()                  
            w_1,reps= func(w_2,k)
            time_line.append(['Enc 1',time.time()-start_1])
            tot_reps+=reps

    if len(w_2)<n:
        start_3 =time.time()
        w_out = enc_3(w_2,n,k)
        time_line.append(['Enc 3',time.time()-start_3])
    else:
        w_out = w_2[:n]
    k_repeats_presents = check_duplicate(w_out,k)
    if k_repeats_presents == False:
        print( "Encoding ended succsesfully ")
    return w_out,tot_reps,time_line


##the decoder control function is in charge of preforming the reverse operations to the encoder in the opposite order
# the last operation of the encoder will be the first to be reversed in the decoder and so on, 
# the last decoder operation will reverse the first encoding operation and will yield the original sequence
# the order, as well as the stoping point, is detemined by the markers and the length of the string. 
def decoder(w_in_dec,k):
    n = len(w_in_dec)
    k = int(np.ceil(2*(1+np.log2(n))))
    if k % 2 != 0:
        k=k+1
   
    w = remove_addons(w_in_dec,n,k)

    while len(w)<n-1:
    # enc 1 shortens the string by one bit per itteration, therefoe if the string is shorter then the original(+prefix marker) there are more enc 1 operation to reverse
        
        if w[0]== 0:  # markes an enc1 operation
            w=dec_1(w,k)

        elif w[0] == 1: # marks an enc2 operation    
            w=dec_2(w,k)

     # length is identical to the original(+prefix marker) there are no more enc 1 operation to reverse 
    while w[0]==1: # markes an enc2 operation (does not changew the length of the string)
        w=dec_2(w,k)
   
    # if the length is n-1 and the first bit is '0' (no enc2 operations to reverse, which would have beeen marked by '1') than we have w ='0'+ original string.
    w_out = w[1:]
    return w_out




### the user interface : menu, operating modes, result analyzer, main 

def display_menu():
    # displaying options for operation 
    print("Operating Modes:")
    print("1. Mode 1- encoding only")
    print("2. Mode 2- decoding only")
    print("3. Mode 3 - full process encoding + decoding")
    print("0. Exit")
    print ("note operating modes 1 and 2 require a digital file where as full process operating mode offfers an option for generating a random sequence of a length of your choise")

## encoding mode: search algorithm: default - dynamic hash based search, input type: file, output type: optional{string, text file }, additional inputs: minimal reads length (see explanation in READ ME)
def encoding_mode():
  
    input_data = input("Enter file path and minimal reads length (if the leangth is not known enter 0):")
    path,gamma=input_data.split()
   
    #translating to binary
    w_in,k = seq_generator(path,2)
    tech_restriction = np.floor((k+1)/2)
   
    # encoding
    (w_out,reps,time_line) = encoder(w_in,0,k)
   
    # verifieng the condition for k (repeatition length), gamma*2 must be larger than k for accurate assembly
    if int(gamma) !=0 :
        if int(gamma) < int(tech_restriction):
            print ("file might be to large for sequencing using chosen sequnecer ")
  
    # saving as file or outputing a string
    print("would you like the encoded string to be saved as a text file ? ")
    out_path =input("if so, enter a file path here. otherwise, enter 0")
    
    if out_path!=0:
        back_to_file(w_out,out_path)
    else:
        return w_out


## decoding mode: input type: file, output type: optional{ string, text file } 

def decoding_mode():
    in_path = input("Enter input file path")
    # translating to binary
    w_in,k = seq_generator(in_path,2)
    #decoding
    w_out = decoder(w_in,k)
    #saving as file or outputing a string
    print("would you like the decoded string to be saved as a file ?")
    out_path =input("if so, enter a file path here. otherwise, enter 0")
    if out_path!=0:
        decoded_file = back_to_file(w_out,out_path)
    else:
        return w_out

## the results analyser calcultes the bit error rate using a simple bitwise comparison of input and output strings as well as the runtime for each section of proccessing 

def results_analyzer(w_in,w_out,start_enc,end_enc,end_dec,reps,time_line):
    
     comp = [abs(w_in[i]-w_out[i])for i in range (len(w_in))]
     tot_errors = sum(comp)+abs(len(w_in)-len(w_out)) # in case the strings dont match in size , this can only happen if the decoding went terribly wrong.
     BER = tot_errors/len(w_in)
   
     runtime_tot = end_dec-start_enc
     runtime_enc = end_enc-start_enc
     runtime_dec = end_dec-end_enc  
     runtime_division = [runtime_enc/runtime_tot,runtime_dec/runtime_tot]
     runtime_precentage = [item*100 for item in runtime_division]

     minutes, seconds = divmod(runtime_enc, 60)
     seconds = "{:.0f}".format(seconds)
     minutes = int(minutes)
     runtime_enc  =f"{minutes} minutes and {seconds} seconds"
     minutes, seconds = divmod(runtime_dec, 60)
     seconds = "{:.0f}".format(seconds)
     minutes = int(minutes)
     runtime_dec =f"{minutes} minutes and {seconds} seconds"
    
     print('The string has been reconstructed with a bit error rate of:', BER,'%')
     print ('Throught encoding',reps,'repetitions were found')
     print(f'Total Encoder runtime = {runtime_enc}')
     print(f'The Encoder runtime is {round(runtime_precentage[0],2)}% of the total runtime')
     print(f'Total Decoder runtime = {runtime_dec}')
     print(f'The Decoder runtime is {round(runtime_precentage[1],2)}% of the total runtime')
     print('The time line in the order of operation is', time_line,'runtimes for each function given in seconds')
     
     operation_runtimes = {}
     total_runtime = 0.0
     for operation, runtime in time_line:
        total_runtime += runtime
        if operation in operation_runtimes:
            operation_runtimes[operation] += runtime
        else:
            operation_runtimes[operation] = runtime
     sizes = [runtime / total_runtime for runtime in operation_runtimes.values()]
     plt.pie(sizes, labels=None, autopct="%1.1f%%", startangle=140)
     plt.legend(operation_runtimes.keys(), title="Operations", loc="best")
     plt.axis("equal")  # Equal aspect ratio ensures that pie is drawn as a circle.
     plt.title("Runtime Distribution")
     plt.show()
    


## test mode: search algorithm: optional{dynamic hash based search, original search algorithm, Booyer moore based search},  input type: optional{file,random sequence (length = n)}, output type: results (= runtimes and bit error rate)


def test_mode():
    # choosingt the type of input
    print("choose one of the following input modes")
    print("1. random sequence")
    print("2. digital file")
    print("0. Exit")
    input_mode = input("Enter your choise here")

    if int(input_mode)== 1:
        print ("choose desired length of sequence")
        input_data = input("enter length")
        n=int(input_data)
        w_in,k= seq_generator(n,1) 

    elif int(input_mode) ==2:
        input_data = input("Enter input file path")
        path=input_data
        w_in,k = seq_generator(path,2)
        
    #choosing a type of search algorithm to use during the encoding process 
    print ("choose search algorithm for encoding")
    print ("1. Dynamic Rabin Karp(defualt option)")
    print("2. Original algorithm ")
    print("3. Boyer-Moore algorithm")
    print("0. Exit")
    alg_num = input("enter your selection here")
    time_line = []
   
    #encoding
    start_enc =time.time()
    (w_enc,reps,time_line_enc) = encoder(w_in,alg_num,k)
    end_enc =time.time()
    k_repeats_in_w_enc = check_duplicate(w_enc,k)
    #decoding
    w_out = decoder(w_enc,k)
    back_to_file (w_out, 'C:\\Users\\ronil\\OneDrive\\IMGed.jpg')


    end_dec = time.time()
    time_line = time_line_enc+[['Decoder',end_dec-end_enc]]
    results_analyzer(w_in,w_out,start_enc,end_enc,end_dec,reps,time_line)

    
def main():
    while True:
        display_menu()
        choice = input("Select an operating mode (0-3): ")

        if choice == '1':
            encoding_mode()
        elif choice == '2':
            decoding_mode()
        elif choice == '3':
            test_mode()
        elif choice == '0':
            print("Exiting the program.")
            break
        else:
            print("Invalid choice. Please select a valid option.")

if __name__ == "__main__":
    main()












    