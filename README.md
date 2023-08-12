# K_Repeat_Free
**encoding for repetition-free sequences**.
In this project I have implemented an encoding and decoding scheme as described in the article “Repeat Free Codes” By O. Elishco, E. Yaakobi, R. Gabrys and M. Medard
This work was motivated by applications in DNA-based storage media.
The reading process of DNA is made up of two parts:
The first part is sequencing, in this process the physical reading takes place, and the order of base pairs is determined, long strand of data can’t be accurately sequenced, so the strand is first divided into shorter overlapping sub-strands, known as ‘reads’, that are than sequenced instead of the complete strand. The output of sequencing is a spectrum of reads.
The second part of the DNA reading process is assembly, in this process we assemble the reads in the spectrum to form the original string using the overlaps between them, in this context the overlap length is denoted by k.
 For this process to yield a unique reconstruction of the original data string we require all k-length overlaps to be unique, i.e., no repeating segments of length k in the strand.
In applications for digital information storage, we can encode the data to a k-repeat free constrained codeword in advance and decode it after it has been sequenced and assembled.




**The Encoding and Decoding algorithms, block diagrams, and a short explanation:**
**The encoder**:
![image](https://github.com/ronilev1n/K_Repeat_Free/assets/141573619/8fecade4-754d-4aab-81e2-3a857d581110)

I chose to divide the encoder into three main parts each in charge of a process specified in the algorithm.
 the original string is first padded on both ends by markers, the prefix marker is a single bit ‘0’, the suffix-marker is ‘1’ followed by   〖log_2〗⁡〖(n)+1〗 ‘0’’s. This padded string is fed to Enc 1 block.
Enc 1: This block’s purpose is to eliminate repetitions and it does so in an iterative manner.
Since information regarding the deleted window (namely location and content) is essential for later decoding, Enc 1 does not identify all the repetitions that appear in the string in one go but rather just one in each iteration. Once the first repetition is found the binary representation of each of the two identical windows indices in the string is added to the string from the right, along with a marker ‘0’, and the first of the two is deleted.
This process essentially creates an entirely new string with different patterns and repetitions, therefore the process repeats operating on the new string as it did on the previous one until it arrives at a string with no repetitive k-Mers. This means that there is no way of determining the number of iterations in advance as it has little to do with the number of repetitions in the input string. This method provides an elegant way of storing the data needed for later recovery in layers of iterative overwrites in a particular order while being capable of handling less straight-forward types of repetitions such as overlapping identical windows, (for example k=4, w=011111111100100, how many repetitions can you find here? The 4-mer ‘1111’ appears 6 times but once the first 4 ‘1’s are deleted w_new= 01111100110 now ‘1111’ only appears twice).
 Note each iteration of Enc 1 shortens the length of the string by one bit as k=2log(n)+2 and each address length is log(n) so in addition to the marker bit, k bits are deleted and k-1 bits are added.
To implement this method an incredibly complicated search is needed, which I will elaborate more on in the next section.
Once the Enc 1 block arrives at a string with no repetitive k-Mers, it passes it forward to be processed by Enc2.
Enc 2: this block is responsible for maintaining the uniqueness of the suffix marker, this is required due to the operation of Enc3, the search here is much more straightforward as the pattern is fixed. 
Once Enc 2 finds the pattern somewhere in the string that is not the end of the string, it deletes the 〖log_2〗⁡〖(n)+1〗  ‘o’’s, which are all the bits of the suffix marker sequence accept the first one, and adds its binary address from the right along with a marker ‘1’ to enable the decoder to identify the two distinct types of eliminations.
Note this process does not change the length of the string.
Enc 3: since each iteration shortens the string by one bit and there is no way of predetermining how many operations Enc 1 will perform, in order to have a uniform length output addition of bits to the end of the string is required. This process of expansion must not create repetitions, the additional bits are meaningless and therefore will be decarded in the first stage of the decoding, the decoder recognizes which part to discard using the suffix marker which is why we have Enc2.

**The decoder:** 

![image](https://github.com/ronilev1n/K_Repeat_Free/assets/141573619/5e9e5e42-24da-4875-b1ad-71a8edbc7182)

I chose to divide the encoder into two main parts, each in charge of reversing the corresponding part of the encoder, Dec 1 is in charge of reversing the operations made by Enc 1 i.e., the elimination of repeated k windows, similarly, Dec 2 is in charge of reversing the operations made by Enc 2 i.e., the elimination of suffix marker prior to the one added at the first stage of encoding.
In order to correctly reverse all encoding operations, the decoder must operate in the opposite order to the encoder, therefore the first operation to reverse will be that of Enc3. 
 As I explained in the section on the encoder the suffix marker’s purpose is to allow the decoder to differentiate between the part of the string that contains relevant data, and the rest which is meaningless and was only added during expansion in Enc 3 to maintain a fixed output length.
Note a string will only enter the expansion process if it had more than k/2 windows found and eliminated in Enc1, strings that have fewer windows detected throughout the encoding process will be longer than n but no longer than n+k/2. By outputting the first n bits the only part being removed are the’0’’s at the end of the suffix marker, which is why we begin the decoding process by checking which of these two cases is relevant for the string at hand, followed by the removal of meaningless bits in according to the relevant case. This stage can be thought of as the reversal of Enc 3.
After we have obtained a string that only contains data relevant to reversing the encoding process, we can begin doing so by using two types of indicators: the length of the string and the first bit.
Since the decoder must reverse the operations of the encoder in the exact opposite order it makes use of these indicators to identify the next operation to reverse (Enc1 or Enc2).
Length indicator: 
The string started out n bits long and had at least one of them removed in the initial stage, there are only two options for the length of the string in the remaining steps of the decoding: the string is either of length n-1 or shorter.
Note that the process that takes place in Enc 1 shortens the string while those that take place in Enc 2 don’t have any effect on the length of the string.
Therefore, having a string that is shorter than the original string + single bit as a prefix marker indicates that some operations that change the length of the string have yet to be reversed. 
 Following this logic having a string that is equal in length to the length of the original string + single bit as a prefix marker, suggests that all the process that changes the length of the string has been reversed, this means that the string is either fully decoded or the only operations left to revers are those of Enc2.
We can differentiate the two cases by checking the first bit. 
First-bit indicator:
During encoding Enc 1 and Enc 2 both used markers of a single bit. Enc1 used ‘0’ and Enc 2 used ‘1’.
the markers were placed at the beginning of the string at the end of each process to allow the decoder to identify its next operation to reverse. 
 Note that the prefix marker is identical to Enc 1 marker (‘0’) since identifying the need to reverse its operations can be done using the length of the string, whereas Enc 2 operations can only be identified by its marker.
Using both indicators gives us four options at the beginning of each operation. 
1.	|w|< n-1 & w [0] = 0  
This means that the decoding is not over + the next step is Dec 1
2.	|w|= n-1 & w [0] = 0
This means that there are no more Enc1 operations to reverse + the next step is not Dec 2, which means the decoding is done.
3.	|w|< n-1 & w [0] =1 
4.	|w| = n-1 & w [0] =1
Both options mean the next step is Dec 2, however, the first option indicates that it’s not the last step, as at least one more length-changing operation (Enc1- Dec1) has yet to be reversed.
Despite the seemingly more complicated process, the decoding process was far less computationally intense than the encoder, as it involves no complicated search.

**OPERATING MODES**
The program supports three types of operations:

**Encoding** only:

![image](https://github.com/ronilev1n/K_Repeat_Free/assets/141573619/4894d75f-b935-49ab-8b64-efb4f48cea8d)

input type: digital file (user provides the path)

Encoders search algorithm: dynamic Rabin Karp hash 

Output: encoded sequence of size length(binary(file)) +2 or a digital file

Additional options: in case the sequencing technology is known in advance users can input the nominal minimal read length given in the sequencers spec to make sure the length of repetition is such that the encoding will provide the desired effect.

**'reads' length restriction**
Generally, in this scheme, the size of repetitions, K, is linked to the size of the input sequence. 
k=2(log_2⁡〖(n)+1〗 )
And for the assembly to be performed with no ambiguities **γ_min≥k+1**
Once the sequence generator converts the file into binary the user interface calls on a function to calculate whether γ_min≥⌊((〖2log〗_2⁡〖(n)+3〗 ))/2⌋ 
Note n in the above calculation includes the two-bit redundancy (n= length(binary(file)) +2), γ_min  is given in base pairs, not bits.
In case the requirement is not met a warning notice is outputted to the interface.

**Decoding** only:
![image](https://github.com/ronilev1n/K_Repeat_Free/assets/141573619/01291b8d-69e2-4e35-92b6-a2108544110a)

Input type: digital file (user provides the path)

output: digital file 

**Full process (Encoding + Decoding + Tests)**:

![image](https://github.com/ronilev1n/K_Repeat_Free/assets/141573619/f788a829-e95c-45b7-9acc-151b5414b798)


Input type: user chooses between file upload (user provides a path) and random sequence (user provides the desired length of sequence)
Encoders search algorithm: users can choose between three of the fastest algorithms I tested for the implementation of the encoder the options are as follows:
	the default dynamic R.K hash
	the original algorithm which I constructed specifically for this type of search (will elaborate further in a later section)
	the Boyer-Moore search algorithm 
 
output: runtime analysis, Bit Error Rate

**search algorithm:**
To improve the runtimes of the most computation-intensive function, Enc 1, responsible for repetition search and elimination, I set out to examine the performance of different search algorithms for pattern detection. Among these, the dynamic hash algorithm proved to be the fastest. 

Additionally, I came up with an original algorithm tailored for repetition search in a constantly changing text, which I have dubbed the “original” search algorithm, I will now explain its method of operation:

in repetition search, we iteratively utilize pattern search algorithms. Each window in the string (except for the last one) is a pattern, the appearance of which we look for in the substring that follows it. 

The algorithm is designed for repetition search, and it uses two methods to reduce the number of both pattern searches and bitwise comparisons in each pattern search. 

The first method is inspired by Boyer Moore’s bad character heuristic, which is not very useful in binary strings. The heuristic I used for this algorithm, which I have dubbed the 'bad suffix rule', allows us to skip whole patterns for comparison in case their suffix does not appear later in the string.
For example, if k = 6 and w = 01011100110011, the first window (= 010111) has the suffix 111, which does not appear later in the string. Therefore, the first, second (101110), third (011100), and fourth (111001) windows do not repeat. Noticing this allows us to skip to the fifth window immediately. This window (110011) does repeat.
 Note that the second window’s suffix 110 does appear later in the string, so we could not have skipped the fifth window.

The second method I used is elimination. Using the sum, as well as the suffix, as indicators, we can eliminate windows from the bitwise comparison in case their indicators are different.
For example, if k = 6, window_i = 111110, and window_j = 111100, in a regular bitwise comparison of windows i and j, we would have to perform five comparisons to identify a difference. However, using a sum indicator where sum_i = 5 and sum_j = 4, we can tell them apart using just one comparison, this is of course at the cost of preprocessing each string .

