#Helpful commands for minion running and data management

MinKNOW writes reads to the local machine. Files will always write to /Library/MinKNOW/data/reads/

##### To check out how your run is going:

Count the number of `fast5` files you have writing to your reads directory.

##### Move reads from local machine over to the cluster

Make a library directory over on the cluster. In house, this should be in `/fh/fast/bedford_t/zika-seq/data`

Once the minion is running and `fast5` files are being written to your machine, you can rsync them over to the cluster. Because the file list is so long, use option `-avz`. The command below removes source files from your local machine once the file has been successfully moved to the destination directory. It will run until the sync has finished, then wait for 60 seconds, then look again to see if new files have come in that should now be rsynced over. If there are directories that have been written in to reads that you _do not_ want synced over, include them in `exclude` option.

    while :; do rsync -avz --remove-source-files --exclude={"pass/","fail/","skip/","tmp/"} /Library/MinKNOW/data/reads/ ablack2@rhino:/fh/fast/bedford_t/zika-seq/data/<library-name> ; sleep 60; done

Sometimes however if you haven't removed all source files from previous runs, you might get files from two different libraries read in to the same directory by MinKNOW. The names will be clearly different, however if there's a lot of files then rsync won't build the argument list if you try to sync specific files with a string matching pattern. In such cases I move the reads to a new local directory to do the sort, and then sync that secondary directory over to the cluster. For example:

`mkdir <temp-dir>`

Do a dry run with `echo`

`find . -type f -name '<string-matching-pattern>' -exec echo mv -v {} <temp-dir> \;`

If everything looks good then run.

`find . -type f -name '<string-matching-pattern>' -exec mv -v {} <temp-dir> \;`

##### Remove empty directories still on your local machine

Depending on how MinKNOW is writing the files locally (and this seems to change with MinKNOW version) you may end up with loads of empty directories once all the `fast5` files have been synced to the cluster. Because removing each one individually would be a pain, you can use these shortcuts.

First, do a dry run and check that things look okay with `find <main-dir> -type d -empty -print`

If that looks all good then delete empty directories with `find <main-dir> -type d -empty -delete`

##### Count how many reads you have, even if they are sorted into a variety of subdirectories

It would be super annoying to pipe ls to wc in every subdirectory that gets made (for instance as MinKNOW writes or post barcode demultiplexing). You can sum the total read counts across all subdirectories using the following command.

`find <parent-dir> -type f -print | wc -l`

You can also specify to count specifically fast5 files with:

`find <parent-dir> -name "*fast5" | wc -l`
