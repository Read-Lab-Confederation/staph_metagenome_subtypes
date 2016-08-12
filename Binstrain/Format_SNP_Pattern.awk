
    {
        n = split( $0, a, "" );
        for( i = 1; i <= n; i++ )
        {
            count[a[i]]++;
            pos[a[i]] = sprintf( "%s%d ", pos[a[i]], i );
        }

        min = "";
        for( x in count )
        {
            if( match( x, "[ACGT]" ) && (min == "" || count[x] < count[min] ) )
                min = x;
        }

        print $0, min, pos[min];

        delete count;
        delete pos;
    }

