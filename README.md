# remphasis
emphasis plugins

## Build the packages

```bash
./build_packages.sh
```

This will build the R source packages in `.\`.

## Install the packages

```bash
./install_packages.sh
```

This will install the R source packages in the default library folder.
You might want to set an alternative library path in `~/.Rprofile`:

```
.libPaths("~/R/mypackages")
```
