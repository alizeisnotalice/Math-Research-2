"""Check the material ownership of each declared multiple contact at every usable moving anchor."""
def install():
    import locate_tomar_five_frontier as curved
    cls=curved.Material
    old_anchor=cls.anchor
    old_forms=cls.formsets
    def anchor(self,a,b):
        result=old_anchor(self,a,b)
        if not hasattr(self,"_event_angles"):self._event_angles={}
        self._event_angles[id(result)]=(a,b)
        return result
    def keys(self,a,b):
        st=self.st;pcs,t,p,rot=st.model(a,b)
        table={}
        for name,end,*_ in pcs:
            if name in st.TIES:
                table[name]=tuple((n,*st.primitive_key(n,end,rot[n])[1:])for n in sorted(st.TIES[name])if n)
        assert set(table)==set(st.TIES)
        return table
    def forms(self,l,h,data,target):
        result=old_forms(self,l,h,data,target)
        if result is not None:
            if not hasattr(self,"_checked_normal_keys"):
                self._checked_normal_keys=keys(self,*self.fixed)
            a,b=self._event_angles[id(data)]
            assert keys(self,a,b)==self._checked_normal_keys,("moving anchor material changes",str(a),str(b))
        return result
    cls.anchor=anchor;cls.formsets=forms
