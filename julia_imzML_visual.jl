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
function SaveBitmapCl( name, pixMap::Array{UInt8,2}, colorTable::Array{UInt32,1} )
    # Get image dimensions
    dim = size( pixMap )
    if length( dim ) != 2
      return 0
    end
    # Normalize pixel values to get a more accurate reading of the image
    min_val = minimum(pixMap)
    max_val = maximum(pixMap)
    pixMap = round.(UInt8, 255 * (pixMap .- min_val) ./ (max_val - min_val))
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
    # Save color table
    write( stream, colorTable )
    if length( colorTable ) < 256
        fixTable = zeros( UInt32, 256 - length( colorTable ) )
        write( stream, fixTable )
    end
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