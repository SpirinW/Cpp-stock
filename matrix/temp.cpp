int main()
{
	unsigned int size=atoi(argv[1]);
	ifstream file;
	file.open("input.txt");
	vector <int> weirdmatrix;
	for(int i = 0, temp; i < size*size; ++i)
	{
		file>>temp;
		weirdmatrix.push_back(temp);
	}
	if(weirdmatrix.size()%size != 0)
	{
		cout<<"Incorrect Matrix\n";
		return 0;
	}
	int m=(weirdmatrix.size())/size;
	vector<vector <int>> a(size, vector <int> (m));
	for(int i = 0; i < size; i++)
	{
		for(int j = 0; j < m; j++)
			a[i][j]=weirdmatrix[i*size + j];
	}
	for(auto row:a)
	{
		for(auto j:row)
			cout<<j<<" ";
		cout<<endl;
	}
	cout<<"Input vector with size "<<size<<":\n";
	vector <int> v;
	for(int i = 0, temp; i < size; ++i)
	{
		cin>>temp;
		v.push_back(temp);
	}
	for(auto i:v)
		cout<<i<<" ";
	cout<<endl;
	vector<int> res;
	for(int i = 0; i < m; ++i)
		res.push_back(0);
	for(int i = 0; i < size; ++i)
	{
		for(int j = 0; j < m; j++)
			res[i]+=v[j]*a[j][i];
	}
	for(auto i : res)
		cout<<i<< " ";
	cout<<endl;
	return 0;
}