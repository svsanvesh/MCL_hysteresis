#!/bin/bash

git status

# Prompt for the commit message
read -p "Enter the commit message: " commit_message

# Perform the Git operations
git pull origin main
git add .
git commit -m "$commit_message"
git push origin main
