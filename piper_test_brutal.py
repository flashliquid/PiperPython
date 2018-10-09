
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import seaborn as sns

# Define plotting functions hsvtorgb and piper
def hsvtorgb( H, S, V ):
    '''
    Converts hsv colorspace to rgb
    INPUT:
        H: [mxn] matrix of hue ( between 0 and 2pi )
        S: [mxn] matrix of saturation ( between 0 and 1 )
        V: [mxn] matrix of value ( between 0 and 1 )
    OUTPUT:
        R: [mxn] matrix of red ( between 0 and 1 )
        G: [mxn] matrix of green ( between 0 and 1 )
        B: [mxn] matrix of blue ( between 0 and 1 )
    '''
    # conversion (from http://en.wikipedia.org/wiki/HSL_and_HSV)
    C = V*S
    #print "C", C

    Hs = H / (np.pi/3)
    #print "Hs", Hs

    X  = C * ( 1 - np.abs( np.mod( Hs, 2.0*np.ones_like(Hs) ) - 1 ) )

    #print X
    N  = np.zeros_like( H )
    # create empty RGB matrices
    R = np.zeros_like( H )
    B = np.zeros_like( H )
    G = np.zeros_like( H )
    # assign values
    h  = np.floor(Hs)
    # h=0
    R[h==0]=C[h==0]
    G[h==0]=X[h==0]
    B[h==0]=N[h==0]
    # h=1
    R[h==1]=X[h==1]
    G[h==1]=C[h==1]
    B[h==1]=N[h==1]
    # h=2
    R[h==2]=N[h==2]
    G[h==2]=C[h==2]
    B[h==2]=X[h==2]
    # h=3
    R[h==3]=N[h==3]
    G[h==3]=X[h==3]
    B[h==3]=C[h==3]
    # h=4
    R[h==4]=X[h==4]
    G[h==4]=N[h==4]
    B[h==4]=C[h==4]
    # h=5
    R[h==5]=C[h==5]
    G[h==5]=N[h==5]
    B[h==5]=X[h==5]
    # match values
    m = V - C
    R = R+m
    G = G+m
    B = B+m
    return( R, G, B )

def piper( dat_piper, plottitle, alphalevel, color ):
    '''
    Create a Piper plot
    INPUT:
        dat_piper: [nx8] matrix, chemical analysis in mg/L
                    order: Ca Mg Na K HCO3 CO3 Cl SO4
        plottitle: string with title of Piper plot
        alphalevel: transparency level of points. If 1, points are opaque
        color: boolean, use background coloring of Piper plot
    OUTPUT:
        figure with piperplot
        dictionary with:
            if color = False:
                cat: [nx3] meq% of cations, order: Ca Mg Na+K
                an:  [nx3] meq% of anions,  order: HCO3+CO3 SO4 Cl
            if color = True:
                cat: [nx3] RGB triple cations
                an:  [nx3] RGB triple anions
                diamond: [nx3] RGB triple central diamond
    '''
    # Basic shape of piper plot
    offset = 0.05
    offsety = offset*np.tan(np.pi/3)
    h = 0.5*np.tan(np.pi/3)
    ltriangle_x = np.array([0, 0.5, 1, 0])
    ltriangle_y = np.array([0, h, 0, 0])
    rtriangle_x = ltriangle_x + 2*offset + 1
    rtriangle_y = ltriangle_y
    diamond_x = np.array([ 0.5, 1, 1.5, 1, 0.5 ]) + offset
    diamond_y = h*( np.array( [ 1, 2, 1, 0, 1 ] ) ) + (offset*np.tan(np.pi/3))
    fig = plt.figure()

    ax = fig.add_subplot(111, aspect='equal', frameon=False, xticks=[], yticks=[])

    ax.plot( ltriangle_x, ltriangle_y, '-k' )
    ax.plot( rtriangle_x, rtriangle_y, '-k' )
    ax.plot( diamond_x, diamond_y, '-k' )

    # labels and title
    sns.set_context("notebook")
    plt.title( plottitle )
    plt.text( -offset,      -offset,  u'$Ca^{2+}$'   , horizontalalignment='left', verticalalignment='center' )
    plt.text( 0.5,          h+offset, u'$Mg^{2+}$'   , horizontalalignment='center', verticalalignment='center' )
    plt.text( 1+offset,     -offset,  u'$Na^+ + K^+$', horizontalalignment='right', verticalalignment='center' )
    plt.text( 1+offset,     -offset,  u'$HCO_3^- + CO_3^{2-}$', horizontalalignment='left', verticalalignment='center' )
    plt.text( 1.5+2*offset, h+offset, u'$SO_4^{2-}$' , horizontalalignment='center', verticalalignment='center' )
    plt.text( 2+3*offset,   -offset,  u'$Cl^-$'      , horizontalalignment='right', verticalalignment='center' )

    # Convert chemistry into plot coordinates
    gram_mol={"Ca":40.078, "Mg":24.305, "Na":22.989768, "K":39.0983, "HCO3":61.01714, "CO3":60.0092, "Cl":35.4527, "SO4":96.0636
}
    valence_abs= {"Ca":2.0, "Mg":2.0, "Na":1.0, "K":1.0, "HCO3":1.0, "CO3":2.0, "Cl":1.0, "SO4":2.0}

    #gmol  = np.array( [40.078, 24.305, 22.989768, 39.0983, 61.01714, 60.0092, 35.4527, 96.0636] )
    #eqmol = np.array( [2., 2., 1., 1., 1., 2., 1., 2. ])


    #n     = dat_piper.shape[0]
    #gram_mol={"Ca":40.078,"Mg":24.305,"Na":22.989768, "K":39.0983,"HCO3":61.01714,"CO3":60.0092,"Cl":35.4527,"SO4":96.0636}
    #valence_abs= {"Ca":2.0, "Mg":2.0, "Na":1.0, "K":1.0, "HCO3":1.0, "CO3":2.0, "Cl":1.0, "SO4":2.0}

    df_dat = df[["Ca", "Mg", "Na", "K", "HCO3", "CO3", "Cl", "SO4"]]
    df_meqL = pd.DataFrame() #initialize an empty dataframe
    for i in gram_mol:
        df_meqL[i]=  (df_dat[i]/gram_mol[i])*valence_abs[i]

    df_meqL = df_meqL[["Ca", "Mg", "Na", "K", "HCO3", "CO3", "Cl", "SO4"]]#reordering(check if this is necessary and or effective)
    sumcat_pd = df_meqL["Ca"] + df_meqL["Mg"] + df_meqL["Na"] + df_meqL["K"] #sum cations
    suman_pd = df_meqL["HCO3"] + df_meqL["CO3"] + df_meqL["Cl"] + df_meqL["SO4"] #sum anions

    '''meqL  = ( dat_piper / gmol ) * eqmol
    sumcat_pd = np.sum( meqL[:,0:4],axis=1)
    suman_pd  = np.sum( meqL[:,4:8],axis=1)
    cat = np.zeros( (n,3) )
    an  = np.zeros( (n,3) )'''
    cat = pd.DataFrame() #initialize an empty dataframe
    an = pd.DataFrame() #initialize an empty dataframe
    cat["Ca"] = df_meqL["Ca"] / sumcat_pd # Ca as proportion
    cat["Mg"] = df_meqL["Mg"] / sumcat_pd # Mg as proportion
    cat["Na + K"] = ( df_meqL["Na"]+df_meqL["K"] ) / sumcat_pd # Na+K as proportion
    an["HCO3 + CO3"]  = (df_meqL["HCO3"]+ df_meqL["CO3"]) / suman_pd # HCO3 + CO3 as proportion
    an["Cl"]  = df_meqL["Cl"] / suman_pd # Cl as proportion
    an["SO4"]  = df_meqL["SO4"] / suman_pd # SO4 as proportion

    # Convert into cartesian coordinates
    cat_x = 0.5*( 2*cat["Na + K"] + cat["Mg"] )
    cat_y = h*cat["Mg"]
    an_x = 1 + 2*offset + 0.5*( 2*an["Cl"] + an["SO4"] )
    an_y  = h*an["SO4"]
    d_x   = an_y/(4*h) + 0.5*an_x - cat_y/(4*h) + 0.5*cat_x
    d_y   = 0.5*an_y + h*an_x + 0.5*cat_y - h*cat_x
    #df_cart = pd.DataFrame()
    #carts_list = carts_list{}
    #carts_list= carts_list.append ("cat_x","cat_y","an_x", "an_y", "d_x", "d_x", "d_y")
    df['cat_x'] = cat_x
    df['cat_y'] = cat_y
    df['an_x'] = an_x
    df['an_y'] = an_y
    df['d_x'] = d_x
    df['d_y'] = d_y
    print (df.head())
    # plot data
    if color is False:
        sns.set()
        sns.scatterplot('cat_x','cat_y', data=df,hue='Location',edgecolor='none', legend=None)#, hue='Location')
        sns.scatterplot('an_x' ,'an_y', data=df, hue='Location',edgecolor='none',legend=None)#, hue='Location', legend=False)
        sns.scatterplot('d_x', 'd_y',data= df, edgecolor='none',hue='Location')#, hue='Location', legend=False)
    else:

        sns.set()
        sns.scatterplot('cat_x','cat_y', data=df,color=".2",edgecolor='k')#, hue='Location')
        sns.scatterplot('an_x' ,'an_y', data=df, color=".2", edgecolor='k')#, hue='Location', legend=False)
        sns.scatterplot('d_x', 'd_y',data= df, color=".2", edgecolor='k')#, hue='Location', legend=False)
    #plt.plot( cat_x, cat_y, '.k',  alpha=alphalevel )
    #plt.plot( an_x,   an_y, '.k', alpha=alphalevel )
    #plt.plot( d_x,     d_y, '.k', alpha=alphalevel )

    # color coding Piper plot
    if color == False:
    # add density color bar if alphalevel < 1
        if alphalevel<1.0:
            ax1  = fig.add_axes([0.75, 0.4, 0.01, 0.2])
            cmap = plt.cm.gray_r
            #cmap = df['Location']
            norm = mpl.colors.Normalize(vmin=0, vmax=1/alphalevel)
            cb1  = mpl.colorbar.ColorbarBase(ax1, cmap=cmap,
                                             norm=norm,
                                             orientation='vertical')
            cb1.set_label('Dot Density')

        return( dict( cat=cat,
                      an=an ) )

        #print dict(cat=cat, an=an)
    else:
        import scipy.interpolate as interpolate
        # create empty grids to interpolate to
        x0 = 0.5
        y0 = x0 * np.tan( np.pi/6 )
        X  = np.reshape( np.repeat( np.linspace(0, 2+2*offset,1000), 1000 ), (1000,1000), 'F' )
        Y  = np.reshape( np.repeat( np.linspace(0, 2*h + offsety ,1000), 1000 ), (1000,1000), 'C' )
        H  = np.nan * np.zeros_like( X )
        S  = np.nan * np.zeros_like( X )
        V  = np.nan * np.ones_like( X )
        A  = np.nan * np.ones_like( X )
        # create masks for cation, anion triangle and upper and lower diamond
        ind_cat = np.logical_or( np.logical_and( X<0.5, Y<2*h*X ),
                                 np.logical_and( X>0.5, Y<(2*h*(1-X) ) ) )
        ind_an  = np.logical_or( np.logical_and( X<1.5+(2*offset), Y<2*h*(X-1-2*offset) ),
                                 np.logical_and( X>1.5+(2*offset), Y<(2*h*(1-(X-1-2*offset)) ) ) )
        ind_ld  = np.logical_and( np.logical_or( np.logical_and( X<1.0+offset, Y > -2*h*X + 2*h*(1 + 2*offset) ),
                                                 np.logical_and( X>1.0+offset, Y >  2*h*X - 2*h ) ),
                                  Y < h+offsety )
        ind_ud  = np.logical_and( np.logical_or( np.logical_and( X<1.0+offset, Y <   2*h*X ),
                                                 np.logical_and( X>1.0+offset, Y <  -2*h*X + 4*h*(1+offset) ) ),
                                  Y > h+offsety )
        ind_d   = np.logical_or( ind_ld==1, ind_ud==1 )

        # Hue: convert x,y to polar coordinates
        # (angle between 0,0 to x0,y0 and x,y to x0,y0)
        H[ ind_cat ] = np.pi + np.arctan2( Y[ind_cat]-y0, X[ind_cat]-x0 )
        H[ ind_cat ] = np.mod( H[ind_cat]-np.pi/6, 2*np.pi )
        H[ ind_an ]  = np.pi + np.arctan2( Y[ind_an]-y0, X[ind_an]- ( x0+1+(2*offset) ) )
        H[ ind_an ]  = np.mod( H[ind_an]-np.pi/6, 2*np.pi )
        H[ ind_d ]   = np.pi + np.arctan2( Y[ind_d]-(h+offsety), X[ind_d]-(1+offset) )
        # Saturation: 1 at edge of triangle, 0 in the centre,
        # Clough Tocher interpolation, square root to reduce central white region
        xy_cat = np.array( [ [0.0,0.0],
                             [ x0,  h],
                             [1.0,0.0],
                             [ x0, y0] ] )
        xy_an  = np.array( [ [ 1+(2*offset)   ,0.0],
                             [ x0+1+(2*offset),  h],
                             [ 2+(2*offset)   ,0.0],
                             [ x0+1+(2*offset), y0] ] )
        xy_d   = np.array( [ [ x0+offset  ,  h+offsety ],
                             [ 1+offset   , 2*h+offsety],
                             [ x0+1+offset,   h+offsety],
                             [ 1+offset   ,     offsety],
                             [ 1+offset   ,   h+offsety] ] )
        z_cat  = np.array( [1.0, 1.0, 1.0, 0.0] )
        z_an   = np.array( [1.0, 1.0, 1.0, 0.0] )
        z_d    = np.array( [1.0, 1.0, 1.0, 1.0, 0.0] )
        s_cat  = interpolate.CloughTocher2DInterpolator( xy_cat, z_cat )
        s_an   = interpolate.CloughTocher2DInterpolator( xy_an , z_an  )
        s_d    = interpolate.CloughTocher2DInterpolator( xy_d  , z_d   )
        S[ind_cat] = s_cat.__call__(X[ind_cat],Y[ind_cat])
        S[ind_an ] = s_an.__call__(X[ind_an] ,Y[ind_an] )
        S[ind_d  ] = s_d.__call__(X[ind_d]  ,Y[ind_d]  )
        # Value: 1 everywhere
        V[ind_cat] = 1.0
        V[ind_an ] = 1.0
        V[ind_d  ] = 1.0
        # Alpha: 1 everywhere
        A[ind_cat] = 1.0
        A[ind_an ] = 1.0
        A[ind_d  ] = 1.0
        # convert HSV to RGB
        R,G,B = hsvtorgb( H, S**0.5, V )
        RGBA = np.dstack( (R,G,B,A) )
        # visualise
        sns.set()
        sns.set_context("talk")
        plt.imshow( RGBA,
                    origin='lower',
                    aspect=1.0,
                    extent=(0,2+2*offset,0,2*h+offsety) )
        # calculate RGB triples for data points
        # hue
        hcat = np.pi + np.arctan2( cat_y-y0, cat_x-x0 )
        hcat = np.mod( hcat-np.pi/6, 2*np.pi )
        han  = np.pi + np.arctan2( an_y-y0, an_x- ( x0+1+(2*offset) ) )
        han  = np.mod( han-np.pi/6, 2*np.pi )
        hd   = np.pi + np.arctan2( d_y-(h+offsety), d_x-(1+offset) )
        # saturation
        scat = s_cat.__call__( cat_x, cat_y )**0.5
        san  = s_an.__call__(   an_x,  an_y )**0.5
        sd   = s_d.__call__(     d_x,   d_y )**0.5
        # value
        v = np.ones_like( hd )
        # rgb
        cat = np.vstack( ( hsvtorgb( hcat, scat, v ) ) ).T
        an  = np.vstack( ( hsvtorgb(  han,  san, v ) ) ).T
        d   = np.vstack( ( hsvtorgb(   hd,   sd, v ) ) ).T
        return( dict( cat     = cat,
                      an      = an,
                      diamond = d ) )


###
# Plot example data
###

# Load data
'''try:
df=pd.read_excel(r path)
except:
df= pd.read_csv(r path)'''


pathx= raw_input("Enter location of data file:")
if len(pathx)<1:
    pathx = '"C:\Users\David\Anaconda2\Progs\GW_data.csv"'
pathx =  pathx.replace( "\\", "/")#avoid escape character issues

pathx = pathx[1:len(pathx)-1] #remove quotes

print pathx

try:
    df=pd.read_excel(pathx)
except:
    df= pd.read_csv(pathx)

title = raw_input("Enter a name for this plot:")
if len(title)<1:
    title = ""
colour = raw_input("colour coded location map? (gives a trippy piper too) y/n: ")
if colour == 'y':
    colour = True
else: colour = False


alpha = 1

#gram_mol={"Ca":40.078,"Mg":24.305,"Na":22.989768, "K":39.0983,"HCO3":61.01714,"CO3":60.0092,"Cl":35.4527,"SO4":96.0636}
#valence_abs= {"Ca":2.0, "Mg":2.0, "Na":1.0, "K":1.0, "HCO3":1.0, "CO3":2.0, "Cl":1.0, "SO4":2.0}

#df_dat = df[["Ca", "Mg", "Na", "K", "HCO3", "CO3", "Cl", "SO4"]]




#dat = np.loadtxt( 'GW_Data.csv',
                  #delimiter=',',
                  #skiprows=1 )
# define extent map
bbox = np.array( [ np.min(df["Latitude"]),
                   np.min(df["Longitude"]),
                   np.max(df["Latitude"]),
                   np.max(df["Longitude"]) ] )
bbox = bbox + np.array( [-0.05*(bbox[2]-bbox[0]),
               -0.05*(bbox[3]-bbox[1]),
                0.05*(bbox[2]-bbox[0]),
                0.05*(bbox[3]-bbox[1]) ] )
'''
*********************
Plot example data
*********************'''



# Piper plot
rgb = piper( df, title, alphalevel=alpha, color=colour )
# Maps (only data points, no background shape files)
print 'rgb', rgb['cat']
'''
*************
In colour
*************
'''
if colour is True:
    sns.set_context("talk")
    fig = plt.figure()
    # cations
    ax1 = fig.add_subplot(131,
                          aspect='equal',
                          facecolor=[0.75,0.75,0.75,1] )
    ax1.scatter( df["Latitude"],
                 df["Longitude"],
                 s=10,
                c=rgb['cat'],
                #c=l,
                 edgecolor='None',
                 zorder = 3)
    plt.title( title +' Cations' )
    plt.xlabel( 'Longitude' )
    plt.ylabel( 'Latitude' )
    plt.grid()
    plt.xlim( bbox[0], bbox[2] )
    plt.ylim( bbox[1], bbox[3] )
    # central diamond
    ax2 = fig.add_subplot(132,
                          aspect='equal',
                          facecolor=[0.75,0.75,0.75,1] )
    ax2.scatter( df["Latitude"],
                 df["Longitude"],
                 s=10,
                 c=rgb['diamond'],
                 #c=n,
                 edgecolor='None',
                 zorder = 3)
    plt.title(title + ' Central Diamond' )
    plt.xlabel( 'Longitude' )
    plt.ylabel( 'Latitude' )
    plt.xlim( bbox[0], bbox[2] )
    plt.ylim( bbox[1], bbox[3] )
    plt.grid()
    # anions
    ax3 = fig.add_subplot(133,
                          aspect='equal',
                          facecolor=[0.75,0.75,0.75,1] )
    ax3.scatter( df["Latitude"],
                 df["Longitude"],
                 s=10,
                 c=rgb['an'],
                 #c=m,
                 edgecolor='None',
                 zorder = 3)
    plt.title( title +' Anions' )
    plt.xlabel( 'Longitude' )
    plt.ylabel( 'Latitude' )
    plt.xlim( bbox[0], bbox[2] )
    plt.ylim( bbox[1], bbox[3] )
    plt.grid()
    sns.set_context("poster",font_scale=5.5)
    plt.show()

#******************
#Not in colour
#******************

else:
    sns.set_context("talk")
    fig = plt.figure()
    # cations
    ax1 = fig.add_subplot(131,
                          aspect='equal',
                          facecolor=[0.75,0.75,0.75,1] )
    ax1.scatter( df["Latitude"],
                 df["Longitude"],
                 s=10,
                #c=rgb['cat'],
                c='.2',
                 edgecolor='None',
                 zorder = 3)
    plt.title( title +' Cations' )
    plt.xlabel( 'Longitude' )
    plt.ylabel( 'Latitude' )
    plt.grid()
    plt.xlim( bbox[0], bbox[2] )
    plt.ylim( bbox[1], bbox[3] )
    # central diamond
    ax2 = fig.add_subplot(132,
                          aspect='equal',
                          facecolor=[0.75,0.75,0.75,1] )
    ax2.scatter( df["Latitude"],
                 df["Longitude"],
                 s=10,
                 #c=rgb['diamond'],
                 c='.2',
                 edgecolor='None',
                 zorder = 3)
    plt.title(title + ' Central Diamond' )
    plt.xlabel( 'Longitude' )
    plt.ylabel( 'Latitude' )
    plt.xlim( bbox[0], bbox[2] )
    plt.ylim( bbox[1], bbox[3] )
    plt.grid()
    # anions
    ax3 = fig.add_subplot(133,
                          aspect='equal',
                          facecolor=[0.75,0.75,0.75,1] )
    ax3.scatter( df["Latitude"],
                 df["Longitude"],
                 s=10,
                 #c=rgb['an'],
                 c='.2',
                 edgecolor='None',
                 zorder = 3)
    plt.title( title +' Anions' )
    plt.xlabel( 'Longitude' )
    plt.ylabel( 'Latitude' )
    plt.xlim( bbox[0], bbox[2] )
    plt.ylim( bbox[1], bbox[3] )
    plt.grid()
    sns.set_context("poster",font_scale=5.5)
    plt.show()
