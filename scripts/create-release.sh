#!/bin/bash
# ARR-MEDIC CYP3A4 Release Creation Script
# Usage: ./scripts/create-release.sh v1.1.0

set -e

VERSION=${1:-"v1.1.0"}
echo "🚀 Creating release for version: $VERSION"

# Validate version format
if [[ ! $VERSION =~ ^v[0-9]+\.[0-9]+\.[0-9]+$ ]]; then
    echo "❌ Invalid version format. Use: v1.1.0"
    exit 1
fi

# Check if CHANGELOG.md exists and contains the version
VERSION_NO_V="${VERSION#v}"
if ! grep -q "\[$VERSION_NO_V\]" CHANGELOG.md; then
    echo "⚠️  Version $VERSION_NO_V not found in CHANGELOG.md"
    echo "📝 Please update CHANGELOG.md first"
    exit 1
fi

echo "✅ Version $VERSION_NO_V found in CHANGELOG.md"

# Create and push git tag
echo "🏷️  Creating git tag: $VERSION"
git tag -a "$VERSION" -m "Release $VERSION

🌐 Multilingual Korean-English CYP3A4 prediction system
🚀 Auto-sync deployment pipeline
📊 Enhanced documentation and examples

See CHANGELOG.md for detailed changes."

echo "📤 Pushing tag to GitHub..."
git push origin "$VERSION"

echo "🎉 Release creation initiated!"
echo "📋 GitHub Actions will automatically:"
echo "   - Extract release notes from CHANGELOG.md"
echo "   - Create GitHub Release"
echo "   - Attach demo assets"
echo ""
echo "🔗 Monitor progress at:"
echo "   https://github.com/Flamehaven/Arr-Medic-CYP3A4/actions"
echo ""
echo "✅ Release will be available at:"
echo "   https://github.com/Flamehaven/Arr-Medic-CYP3A4/releases"