#https://gist.github.com/jennybc/362f52446fe1ebc4c49f

#Suppress warnings here since if running locally, the /active path does not exist.
RPROJ <- suppressWarnings(list(PROJHOME = normalizePath(getwd()),
                               ACTIVE=normalizePath("/active"),
                               ARCHIVE=normalizePath("/archive")))

host <- Sys.info()[["nodename"]]
local_hostname <- "MLVKQ0PF47H0"
if (length(dir("~/active/")) > 0 & host == local_hostname){
  #Note: on my local machine, these only work when connected to VPN and Then `smbf_mount` mounted the network drives for /active and /archive
  RPROJ[c("ACTIVE","ARCHIVE")] <- paste0("~",RPROJ[c("ACTIVE","ARCHIVE")])

}

attach(RPROJ)
rm(RPROJ,host,local_hostname)
