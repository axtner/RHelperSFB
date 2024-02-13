test=readLines("/home/bioadmin/Projects/BiDoup_2022/221114_M01108_0115/SampleSheet.csv")
test[1]
test[test=="[Reads]"]
length(test)
for(n in 1:length(test)){
  if(test[n] == "[Reads]"){
    print(n)
    line= n
  }
}

read_no = test[line+1]
