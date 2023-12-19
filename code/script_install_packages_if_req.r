packs <- c("XML", "devtools", "RCurl", "fakePackage", "SPSSemulate")
success <- suppressWarnings(sapply(packs, require, character.only = TRUE))
success == F
names(success[success == F])
install.packages(names(success)[!success])
sapply(names(success)[!success], require, character.only = TRUE)
names(success)
