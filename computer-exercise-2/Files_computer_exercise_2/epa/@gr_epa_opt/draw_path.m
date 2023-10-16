function draw_path( gp )
hold( 'on' );

h = plot( gp.path( :, 2 ), gp.path( :, 3 ), ...
	  sprintf( '%sx-', gp.path_color ) );
      
set( h, 'HitTest', 'off', 'HandleVisibility', 'on', ...
	'LineWidth', 2 );
hold( 'off' );