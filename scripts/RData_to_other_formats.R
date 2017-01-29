#
# Transfrom RData File to Other Formats
#
# Based on a script by Maren Buttner.
#
# Copyright (c) 2016 F. Alexander Wolf (http://falexwolf.de).
#


# specify the RData file name
datafilename = "Paul_Cell_MARSseq_GSE72857.RData"
# specify the name of the outfile without extension
outbasename = "paul15"

# choose format: must be one of "csv", "xlsx" (stores multiple files, single csv
# and or xlsx files per data object) or "xlsx", "onehdf5" (stores one file for
# all data frames)
#format = "csv"
format = "onehdf5"
#format = "onexlsx"

if (format == "onehdf5")
{
    library(rhdf5)
    outfilename <- paste(outbasename,'.h5',sep="") 
    h5createFile(outfilename)
}

if (format == "onexlsx")
{
    # the xlsx package needs java
    # library(xlsx)
    # the openxlsx package doesn't need java, and provides some additional
    # functionality, like creating workbooks to write to a common single excel file
    library(openxlsx)
    wb <- createWorkbook()
    outfilename <- paste(outbasename,'.xlsx',sep="") 
}

objlist = load(datafilename)

for(objname in objlist)
{
    print(paste("-----",objname,"-----"))
    # as data frame
    if(format!="onehdf5")
        obj <- as.data.frame(get(objname))
    # as matrix, h5write seems to have problems with some data frames
    if(format=="onehdf5")
        obj <- get(objname)
    # print properties of the object
    # print(head(obj,n=3))
    # print(class(obj))
    # print(dim(obj))
    # print(head(rownames(obj)))
    # print(!is.null(rownames(obj)))
    # print(head(colnames(obj)))    
    # print(!is.null(colnames(obj)))

    if(format=="csv")
        write.csv(obj,file=paste(outbasename,"_",objname,".csv",sep="")) 

    if(format=="onehdf5")
    {

        # write rownames and colnames separately
        if(is.matrix(obj))
        {
           rn <- rownames(obj)
           if(!is.null(rn))
               h5write(rn, file = outfilename, 
                       paste(objname,"_rownames",sep=""))
           cn <- colnames(obj)
           if(!is.null(cn))
               h5write(cn, file = outfilename, 
                       paste(objname,"_colnames",sep=""))
        }

        # treat factor vectors separately
        # otherwise this will give errors
        if(is.factor(obj))
        {
            h5write(as.integer(obj), file=outfilename, paste(objname,"_codes",sep=""))
            h5write(levels(obj), file=outfilename, paste(objname,"_strings",sep=""))    
        }
        # write other matrices and vectors here
        else
            # here we assume we will read with C-style program
            # where data looks transposed as compared to the Fortran-style R
            # therefore, we therefore transpose the data to match 
            # the C convention
            h5write(t(obj), file = outfilename, objname)

    }

    if(format=="xlsx")
        write.xlsx(df, file=paste(outbasename,"_",objname,".xlsx",sep="")) 
    
    if(format=="onexlsx")
    {
         addWorksheet(wb, objname)
         writeData(wb, objname, x = obj)
    }

}

if(format=="onehdf5")
{
    # see properties of the file
    print(h5ls(outfilename))
    H5close()
}

if(format=="onexlsx")
    saveWorkbook(wb, file = outfilename, overwrite=TRUE)


