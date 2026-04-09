#!/usr/bin/env bash
set -euo pipefail

echo "Restoring archived test data and resources..."
tar -xzf test/compressed_resources/test_data.tar.gz
tar -xzf test/compressed_resources/test_resources.tar.gz

echo "Fetching LOLACore resources..."

wget http://big.databio.org/regiondb/LOLACoreCaches_180412.tgz
tar -xzvf LOLACoreCaches_180412.tgz
mv nm/t1/resources/regions/LOLACore/ test/resources/
rm -rf nm
rm LOLACoreCaches_180412.tgz
