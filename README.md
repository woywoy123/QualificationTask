# Qualification Task Attack Plan 

### Next Steps:
1. ~~Verify that RooFit is correctly implemented by building a PDF model and then using the same PDF as data with some applied scaling. For verification to be successful, the fit needs to return this exact scaling number. This proves that the fit is successfully implemented and we can be confident this wont cause any issues.~~ 
2. ~~use the monte carlo and create a sum of the pure n-track templates histogram. Then perform a subtraction of each n-pure to the created histogram and output the result with a fit. It is expected that which ever n-pure was subtracted from the histogram their contribution as measured by the fit would be 0. This verifies two things, the subtraction of the histogram works and the fit was able to resolve the change.~~ 
3. ~~Repeat step 1 but this time with using the monte carlo data to create n-track samples that can have arbitrary truth. Perform the fit and compare the number of entries with the prediction of the fit. The criteria for successful implementation is that the fit is able to exactly predict each n-track classification for the given n-track histogram distribution.~~ (Sort of worked well to approximately a few percent) 
4. Repeat step 2 with 4-tracks histogram and assume this sample to not contain any other cross track contamination. You subtract this "4 track template" from the 3 track histogram. Important to remember is to normalize the entries in the 4-track template to the 3-track template in the range 10 -> 20. Then a fit is performed on the resultant distribution. It is expected that most of the 4-track entries are removed from the 3-track template as predicted by the fit. The output of the fit needs to be compared to the monte carlo as to get closure testing. (Close to being completely implemented. Need to do more tests) 
5. Repeat step 4 for 1-track and 2-track. (Attempted but seems to be troubling. Need to find a proper domain to do the subtraction on) 



