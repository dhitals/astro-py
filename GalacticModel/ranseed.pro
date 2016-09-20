function ranseed
; Returns a "random" long integer.
;
; Pick a random integer between 0001 and 10000 by looking at the
; system clock.
; Useful for initializing the random number generator, since IDL
; doesn't pick good values when seed in undefined.
;
; Obtained from Phil Everson, 12/11/01
;
; Phil says:
;
; Hi Eric,
; Andrea and I have both been frustrated by IDL's seed function. As you 
; noticed, it seems to advance the seed by only one increment for each 
; call (even if you have it generate 1000's of random values). I wrote 
; this function as a quick fix.
; The function ranseed calls systime to get the number of seconds since 
; some long ago moment. That number doesn't change very fast, so I 
; reorder the digits, beginning with 10th's of a second, followed by 
; seconds, 10's of seconds and 100's of seconds to get a 4-digit 
; integer. Then I run this through a linear congruential generator 20 
; times to mix through the possible long-integers. The result is a 
; long-integer chosen from one of 10,000 widely dispersed values, and I 
; use this as my starting seed. Within another program calling ranu, I 
; first give the command  seed=1L*ranseed()  to intialize the seed 
; value. John Bocio just gave me a program that does something sneakier 
; to get a starting seed, but I don't yet have a version for IDL.
;
; Phil

secs = systime(1)
secs = secs - floor(secs/1000.0)*1000
startseed = floor(secs*10.0)
; startseed is the 4-digit number made up of the 100's, 10's,
; 1's and tenths digits returned by systime.
d1 = startseed/1000
d2=(startseed-1000*d1)/100
d3=(startseed-1000*d1 - 100*d2)/10
d4=(startseed-1000*d1 - 100*d2 - 10*d3)
newseed=d4*1000L+d3*100L + d2*10L + d1 + 1L
; newseed is constructed by reversing the order of the
; digits in startseed and adding 1 (to prevent seed=0).
;
; Now run newseed through the minimal standard linear
; congruential generator 20 times, to shuffle through
; long integer values.
;
a = 16807.0d
m = 2147483647.0d
k=1.0*newseed
rep = 20
for i = 1, rep do begin
    c=k*a
    d=floor(c/m)*1.0d
    k=c-d*m
endfor
seed = k*1L
return, seed
end
