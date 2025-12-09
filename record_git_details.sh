# Get current commit hash
COMMIT_HASH=$(git rev-parse HEAD)
echo "Running with commit: $COMMIT_HASH" > version.txt

# Get branch name
BRANCH_NAME=$(git rev-parse --abbrev-ref HEAD)
echo "On branch: $BRANCH_NAME" >> version.txt

# Get status (modified files, etc.)
git status --short >> version.txt

git diff >> version.txt
