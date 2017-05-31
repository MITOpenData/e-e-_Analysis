void exercise(){


 TVector3 myvector ( 0, 0, 1);
 Double_t carray[3];
 myvector.GetXYZ( carray );
 
cout << "x = " << carray[0] << endl;
cout << "y = " << carray[1] << endl;
cout << "z = " << carray[2] << endl;

TVector3 direction (1,0,0);
direction.GetXYZ( carray );

cout << "x = " << carray[0] << endl;
cout << "y = " << carray[1] << endl;
cout << "z = " << carray[2] << endl;


myvector.RotateUz(direction);
myvector.GetXYZ( carray );
 
cout << "x = " << carray[0] << endl;
cout << "y = " << carray[1] << endl;
cout << "z = " << carray[2] << endl;







}
