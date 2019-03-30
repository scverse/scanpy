convert "$1" \
    \( -clone 0 -background none -gravity center -resize  16x -extent  16x16 \) \
    \( -clone 0 -background none -gravity center -resize  32x -extent  32x32 \) \
    \( -clone 0 -background none -gravity center -resize  48x -extent  48x48 \) \
    \( -clone 0 -background none -gravity center -resize  64x -extent  64x64 \) \
    \( -clone 0 -background none -gravity center -resize 128x -extent 128x128 \) \
    \( -clone 0 -background none -gravity center -resize 256x -extent 256x256 \) \
    -delete 0 favicon.ico