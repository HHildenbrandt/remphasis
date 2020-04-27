locate_plugin = function(PlugIn) {
  path = getLoadedDLLs()[[PlugIn]][["path"]]
  if (!is.character(path)) {
    msg = paste(c("Can't locate plugin. Did you load library '", PlugIn, "'?"), collapse = '')
    warning(msg)
    return(PlugIn)  # could be absolute path
  }
  return(path)
}
