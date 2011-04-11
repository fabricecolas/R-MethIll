summary.MethIll <-
function(object, ...){
   rd <- print(object)
   cat('\nSamples with Z.K.IAC>-2\n')
   print(rd[rd$Z.K.IAC>-2,])
   cat('\nSamples with Z.K.IAC<=-2 and Z.K.IAC>-3\n')
   print(rd[rd$Z.K.IAC<= -2 & rd$Z.K.IAC> -3,])
   cat('\nSamples with Z.K.IAC<=-3\n')
   print(rd[rd$Z.K.IAC<= -3,])
}

