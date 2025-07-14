echo "Running run.sh script"

python3 src/main.py --no-heatmap --quiet

echo "Exported to output/clusters_export.json"

echo "filtering (truncation)"

python3 cluster_filter.py output/clusters_export.json output/intermediate.json --mode truncate --max-length 80

echo "filtering (removal)"

python3 cluster_filter.py output/intermediate.json output/filtered.json --mode remove --pretty

echo "filtering (removing features)"

python3 cluster_filter.py output/filtered.json output/filtered.json --mode remove --fields "features" --pretty

echo "Fully filtered!"
