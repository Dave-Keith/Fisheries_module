stretch.beta.val <- function(mn,sd,minb,maxb,fx)
{
# Now this is a stretched beta distribution, this one is used for
# fertility and stuff when you need a beta distribution that isn't
# bounded between 0 and 1.  Better than the log.normal!!

if(sd == 0) bb <- mn # i.e. no variaion therefor the value is just = to the mean!!!
   # Convert the beta parameters to corresponding ones fro a (0,1) beta distribution
 if(sd>0)
{
    mn.beta <- (mn-minb) / (maxb-minb)
    sd.beta <- sd/(maxb-minb)
    # now make sure that the parameter combos are actually possible
if (sd.beta <= ((mn.beta*(1-mn.beta))^0.5))
{
    b.value <- beta.val(mn.beta,sd.beta,fx) # This is my beta function "beta_dist_rand_num.R"
    bb <- b.value * (maxb-minb) + minb # This is the conversion to a stretched beta
} # end the second if

# Now I want to revisit this function and make it so I adjust the sd to the maximum possible sd
# that way we can still use this function when actual sd is too high for the function...
 
 if (sd.beta > ((mn.beta*(1-mn.beta))^0.5))
 {
     sd.beta <- 0.95*((mn.beta*(1-mn.beta))^ 0.5)*(maxb-minb)
     b.value <- beta.val(mn.beta,sd.beta,fx) # This is my beta function "beta_dist_rand_num.R"
     bb <- b.value * (maxb-minb) + minb # This is the conversion to a stretched beta
     message(paste("WARNING!!! your sd was too high"))
     message(paste("I have adjusted it to 95% of max allowable"))
     max.sd <- ((mn.beta*(1-mn.beta))^ 0.5)*(maxb-minb)
     message(paste("please check your output to ensure this give reasonable values"))
     message(paste("Have a nice day!"))
     } # end the second if
 
#else
#{
#    message(paste("Well jerk, looks like you fucked this up big time"))
#    message(paste("do you even have a clue what you are doing, the sd"))
#    message(paste(" is too high for the mean with a vital rate that has"))
#    message(paste("a mean, sd, and min and max values of"))
#    message(paste(mn, sd, minb, maxb))
#    message(paste("Like don't you know that the max possible sd in this case is"))
#    max.sd <- ((mn.beta*(1-mn.beta))^ 0.5)*(maxb-minb)
#    message(paste(max.sd))
#    message(paste("I'm so outa here..."))
#    break
#}) # this is the end of the second else

 } # end the sd > 0 bit...

return(bb)
} # end the function
