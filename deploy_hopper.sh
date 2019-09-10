DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

LATEST_CONTAINER_BUILD="$( ls -t $DIR/containers/container_*.sif |head -n1)"
CONTAINER_BASENAME=${LATEST_CONTAINER_BUILD##*/}
PIPELINE_DEST="/media/hopper/pipelines/exome"
CONTAINER_DEST=/media/hopper/resources/containers/$CONTAINER_BASENAME


# Deploy container if it isn't already deployed
if test -f "$CONTAINER_DEST"; then
    echo "Latest container already deployed, skipping!"
else
    echo "Deploying container"
    cp $LATEST_CONTAINER_BUILD $CONTAINER_DEST
    # TODO: Replace "active" container symlink on hopper!
fi


# Copy pipeline script
cp $DIR/exome.nf $PIPELINE_DEST

# Copy configuration file
cp $DIR/configs/nextflow.hopper.config $PIPELINE_DEST/nextflow.config

# Copy other files
cp $DIR/test_run_exome.sh $PIPELINE_DEST




