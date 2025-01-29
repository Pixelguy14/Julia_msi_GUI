# IntQuantCl is originally a function of julia mzMl imzML with the adition 
# of altering the scale of the colors according to colorlevel.
function IntQuantCl( slice , colorLevel)
    # Compute scale factor for amplitude discretization
    lower  = minimum( slice )
    scale  = colorLevel / maximum( slice )
    dim    = size( slice )
    image  = zeros( UInt8, dim[1], dim[2] )
    for i in 1:length( slice )
        image[i] = convert( UInt8, floor( slice[i] * scale + 0.5 ) )
    end
    return image
end

# SaveBitmap is also originally a function of the mzML imzML library in julia,
# This function dinamically adjust the color palete adjusting to the ammount of colors
# available in pixmap
function SaveBitmapCl( name,
    pixMap::Array{UInt8,2},
    colorTable::Array{UInt32,1} )
    # Get image dimensions
    dim = size( pixMap )
    if length( dim ) != 2
      return 0
    end
    # Compute row padding
    padding = ( 4 - dim[1] & 0x3 ) & 0x3
    # Compute file dimensions. Header = 14 + 40 + ( 256 * 4 ) = 1078
    offset   = 1078
    imgBytes = dim[2] * ( dim[1] + padding )
    # Create file
    stream = open( name, "w" )
    # Save file header
    write( stream, UInt16( 0x4D42 ) )
    write( stream, UInt32[ offset + imgBytes, 0 , offset ] )
    # Save info header
    write( stream, UInt32[ 40, dim[1], dim[2], 0x80001, 0 ] )
    write( stream, UInt32[ imgBytes, 0, 0, 256, 0 ] )
    # Get unique colors displayed in pixMap
    unique_colors = unique(vec(pixMap))
    n_colors = length(unique_colors)
    # Color levels can't surpass 256
    levels = min(n_colors, 256)
    extendedColorTable = zeros(UInt32, 256)
    for i in 1:256
      idx = ceil(Int, i * levels / 256)
      extendedColorTable[i] = colorTable[idx]
    end
    write(stream, extendedColorTable)
    # Save image pixels
    if padding == 0
        for i = 1:dim[2]
            write( stream, pixMap[:,i] )
        end
    else
        zeroPad = zeros( UInt8, padding )
        for i in 1:dim[2]
            write( stream, pixMap[:,i] )
            write( stream, zeroPad )
        end
    end
    # Close file
    close( stream )
end