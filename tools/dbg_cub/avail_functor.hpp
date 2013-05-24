
template<
  typename Factory_,
  typename Functor_,
  bool variadic_ = (Factory_::variadic != 0)>
class AvailFunctorHelper;

template<
  typename Factory_,
  typename Functor_>
class AvailFunctorHelper<Factory_, Functor_, false>
{
private:
  Functor_& _functor;

public:
  explicit AvailFunctorHelper(Functor_& functor) :
    _functor(functor)
  {
    _functor.add_name(Factory_::name());
  }

  void alias(const String& name)
  {
    _functor.add_alias(name, Factory_::name());
  }
};

template<
  typename Factory_,
  typename Functor_>
class AvailFunctorHelper<Factory_, Functor_, true>
{
private:
  Functor_& _functor;

public:
  explicit AvailFunctorHelper(Functor_& functor) :
    _functor(functor)
  {
    _functor.add_name(Factory_::name() + ":<"
      + stringify(int(Factory_::min_points)) + "-"
      + stringify(int(Factory_::max_points)) + ">");
  }

  void alias(const String& name, Index num_points)
  {
    _functor.add_alias(name, Factory_::name() + ":" + stringify(num_points));
  }
};

class AvailFunctor
{
private:
  StringSet& _names;
  StringMap& _aliases;

public:
  explicit AvailFunctor(StringSet& names, StringMap& aliases) :
    _names(names),
    _aliases(aliases)
  {
  }

  template<typename Factory_>
  void factory()
  {
    AvailFunctorHelper<Factory_, AvailFunctor> functor(*this);
    Factory_::alias(functor);
  }

  void add_name(const String& name)
  {
    _names.insert(name);
  }

  void add_alias(const String& alias, const String& name)
  {
    _aliases.insert(std::make_pair(alias, name));
  }
};
